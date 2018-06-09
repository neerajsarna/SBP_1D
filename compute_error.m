function [] = compute_error(M_values,n,n_ref, ...
                                                    t_end,filename,foldername)
    

    num_samples = length(M_values);
    
    result = cell(num_samples,1);
    
% loop over all the systems 
for i = 1:num_samples
    output_filename = strcat(foldername,'/',filename,...
                    '_tend_',num2str(t_end),'_points_',num2str(n),'_neqn_');
    output_filename = strcat(output_filename,num2str(M_values(i)),'.txt');
    
    result{i} = dlmread(output_filename,'\t');
    delta_x = result{i}(1,2) - result{i}(1,1);
    
    % we remove the contribution the rows containing the spatial locations
    result{i} = result{i}(2:end,:)';
    
end    
    % we need the P matrix for computing the integrals
    [~,P,~] = sbp_traditional_2(delta_x,n);
    
    
%% read the reference solutions
    output_filename = strcat(foldername,'/',filename,'_tend_',...
                        num2str(t_end),'_points_',num2str(n),'_neqn_');
    output_filename = strcat(output_filename,num2str(n_ref),'.txt');
    
    result_ref = dlmread(output_filename,'\t');
    
    result_ref = result_ref(2:end,:)';
    errorTot = zeros(num_samples,1);
    
%% computing the error 
disp('error computation ...');
for i = 1:num_samples
    
    % num rows we need to extract from the higher order moment solution
    num_cols = size(result{i},2);
    
    error = [(result_ref(:,1:num_cols)-result{i}), result_ref(:,num_cols+1:end)];
    
    % compute the error. The quantity dot(error,P*error,1) gives us ||.||^2_{L^2(\Omega;\mbb R^n)}
    errorTot(i) = sqrt(sum(dot(error,P*error,1)));
    
end

   

   %% compute the expected rate of convergence
   
    disp('regularity of reference ...');
   [expected_rate] = expected_convg_rate(n_ref, ...
                                    t_end,n,foldername,filename);

   M_values_O = find(mod(M_values,2) ~= 0);
   M_values_E = find(mod(M_values,2) == 0);
   
   [P_error_O,y_error_O] = polyfit_linear(log(M_values(M_values_O)), ....
                                    log(errorTot(M_values_O))');
    convg_rate_O = abs(P_error_O(1));
    
    [P_error_E,y_error_E] = polyfit_linear(log(M_values(M_values_E)), ...
                                        log(errorTot(M_values_E))');
    convg_rate_E = abs(P_error_E(1));
   %% plotting 
   figure(5);
   loglog(M_values,errorTot,'-o', ...
          M_values(M_values_O),exp(y_error_O),'k-', ...
          M_values(M_values_E),exp(y_error_E),'k-', ...
          'MarkerSize',3);
   xlabel('M','FontSize',10);
   ylabel('||E_M||', 'FontSize',10);
   xlim([M_values(1) M_values(end)]);
   legend('error','linear fit');
    xticks(10:30:100);
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
    grid on;
    
    disp('OBTAINED RATE Odd: ');
    disp(convg_rate_O);
    
    disp('OBTAINED RATE Even: ');
    disp(convg_rate_E);
    
end

function [P,yfit] = polyfit_linear(x,y)
P = polyfit(x,y,1);
yfit = P(1) * x + P(2);
end
    

