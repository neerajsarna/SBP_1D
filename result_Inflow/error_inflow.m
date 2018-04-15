clear all;

% we compute the error from 
% number of intervals
n = 300;

% spacing between two grid points
delta_x = 1/n;

% end time
t_end = 0.3;

% we need the P matrix for computing the error
[D,P,InvP] = sbp_traditional_2(delta_x,n);

% 4-150 and 1 for the reference i.e. 200
num_samples = 97 + 1;

result = cell(num_samples,1);
counter = 1;

M_values = 4:100;

% loop over all the systems 
for i = M_values
    output_filename = strcat('inflow_tend_',num2str(t_end),'_points_',num2str(n),'_neqn_');
    output_filename = strcat(output_filename,num2str(i),'.txt');
    
    result{counter} = dlmread(output_filename,'\t');
    
    counter = counter + 1;
end

    output_filename = strcat('inflow_tend_',num2str(t_end),'_points_',num2str(n),'_neqn_');
    output_filename = strcat(output_filename,num2str(200),'.txt');
    
    result{counter} = dlmread(output_filename,'\t');
    

errorTot = [];

%% computing the error 
for i = 1:num_samples-1
    
    % num rows we need to extract from the higher order moment solution
    num_rows = size(result{i},1);
    
    error = [(result{num_samples}(1:num_rows,:)-result{i}); result{num_samples}(num_rows+1:end,:)];
    
    % we need to transpose to allow multiplication with P
    error = error';

    % compute the error. The quantity dot(error,P*error,1) gives us ||.||^2_{L^2(\Omega;\mbb R^n)}
    errorTot = [errorTot sqrt(sum(dot(error,P*error,1)))];
    
end
%% decay of moments in f
odd_Moments = 2:2:200;
even_Moments = 1:2:200;

Odd_Moments = result{num_samples}(odd_Moments+1,:)';
Even_Moments = result{num_samples}(even_Moments+1,:)';

norm_Odd_f = sqrt(dot(Odd_Moments,P * Odd_Moments,1));
norm_Even_f = sqrt(dot(Even_Moments,P * Even_Moments,1));

norm_Odd_dx_f = sqrt(dot(D * Odd_Moments,P * D * Odd_Moments,1));
norm_Even_dx_f = sqrt(dot(D * Even_Moments,P * D * Even_Moments,1));



%% plotting 

figure(1)
loglog(M_values,errorTot,'-o','MarkerSize',3);
xlabel('M','FontSize',10);
ylabel('E_M', 'FontSize',10);
xlim([3 110]);

xticks(10:30:100);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16);
grid on;

figure(2)
loglog(odd_Moments,norm_Odd_f,'-o',even_Moments,norm_Even_f,'-o','MarkerSize',3);
xlabel('m','FontSize',10);
ylabel('N_m', 'FontSize',10);

xticks(10:30:100);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16);
grid on;

figure(3)
loglog(odd_Moments,norm_Odd_dx_f,'-o',even_Moments,norm_Even_dx_f,'-o','MarkerSize',3);
xlabel('m','FontSize',10);
ylabel('N_m^{x}', 'FontSize',10);

xticks(10:30:100);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16);
grid on;

%% computing the convergence rates
disp('rates for error');
disp(convg_rate(M_values(1:25),errorTot(1:25)));
disp(convg_rate(M_values(26:50),errorTot(26:50)));
disp(convg_rate(M_values(51:num_samples-1),errorTot(51:num_samples-1)));

disp('rates for Odd f');
disp(convg_rate(odd_Moments(1:25),norm_Odd_f(1:25)));
disp(convg_rate(odd_Moments(26:50),norm_Odd_f(26:50)));
disp(convg_rate(odd_Moments(51:100),norm_Odd_f(51:100)));

disp('rates for Even f');
disp(convg_rate(even_Moments(1:25),norm_Even_f(1:25)));
disp(convg_rate(even_Moments(26:50),norm_Even_f(26:50)));
disp(convg_rate(even_Moments(51:100),norm_Even_f(51:100)));

function f = convg_rate(x_values,y_values)
f = log(y_values(end)/y_values(1))/log(x_values(end)/x_values(1));
end