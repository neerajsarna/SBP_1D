function [convg_rate_Odd,convg_rate_Even,expected_rate] = compute_error2D(M_values,n,n_ref,n_ref2,...
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
   [expected_rate,loc_truncate] = expected_convg_rate2D(n_ref,n_ref2,t_end,n, ...
                                                        foldername,filename);
   
% 
   disp('truncated at');
   disp(loc_truncate); 
   
   % order for even values of M
   [P_error_Even,y_error_Even] = polyfit_linear(log(M_values(mod(M_values,2)==0)), ...
                                      log(errorTot(mod(M_values,2)==0))');
    
   % order for odd values of M
   [P_error_Odd,y_error_Odd] = polyfit_linear(log(M_values(mod(M_values,2)~= 0)), ...
                                      log(errorTot(mod(M_values,2)~= 0))');
                                  
   convg_rate_Odd = abs(P_error_Odd(1));
   convg_rate_Even = abs(P_error_Even(1));
    
   %% plotting 
   figure(5);
   loglog(M_values,errorTot,'-o', ...
         M_values(mod(M_values,2)==0),exp(y_error_Even),'k-*', ...
         M_values(mod(M_values,2)~=0),exp(y_error_Odd),'r-*', ...
         'MarkerSize',3);
        
   xlabel('M','FontSize',10);
   ylabel('||E_M||','FontSize',10);
   xlim([M_values(1) M_values(end)]);
   legend('error','linear fit Even','linear fit Odd','Location','best');
    xticks(3:4:27);
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
    grid on;
    
    disp('OBTAINED RATE Odd M: ');
    disp(convg_rate_Odd);
    
    disp('OBTAINED RATE Even M: ');
    disp(convg_rate_Even);
    
    disp(' Error in prediction Odd M');
    disp((convg_rate_Odd-expected_rate));
   
    disp('Error in prediction Even M');
    disp((convg_rate_Even-expected_rate));
end

function [P,yfit] = polyfit_linear(x,y)
P = polyfit(x,y,1);
yfit = P(1) * x + P(2);
end
    

