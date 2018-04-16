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

% 4-150 and 2 for the reference i.e. 200 and 199
num_samples = 97 + 2;

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
    
    result{counter+1} = dlmread(output_filename,'\t');
    

    output_filename = strcat('inflow_tend_',num2str(t_end),'_points_',num2str(n),'_neqn_');
    output_filename = strcat(output_filename,num2str(199),'.txt');
    
    result{counter} = dlmread(output_filename,'\t');
    
errorTot = [];

%% computing the error 
for i = 1:num_samples-2
    
    % num rows we need to extract from the higher order moment solution
    num_rows = size(result{i},1);
    
    error = [(result{num_samples}(1:num_rows,:)-result{i}); result{num_samples}(num_rows+1:end,:)];
    
    % we need to transpose to allow multiplication with P
    error = error';

    % compute the error. The quantity dot(error,P*error,1) gives us ||.||^2_{L^2(\Omega;\mbb R^n)}
    errorTot = [errorTot sqrt(sum(dot(error,P*error,1)))];
    
end
%% decay of moments in f
[norm_Odd_f_200,norm_Even_f_200] = norm_OddEven(result{num_samples},200,P);
[norm_Odd_f_199,norm_Even_f_199] = norm_OddEven(result{num_samples-1},199,P);

[norm_Odd_f_dx_200,norm_Even_f_dx_200] = norm_dx_OddEven(result{num_samples},200,P,D);
[norm_Odd_f_dx_199,norm_Even_f_dx_199] = norm_dx_OddEven(result{num_samples-1},199,P,D);

%% plotting 

figure(1);
loglog(M_values,errorTot,'-o','MarkerSize',3);
xlabel('M','FontSize',10);
ylabel('E_M', 'FontSize',10);
xlim([3 110]);

xticks(10:30:100);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16);
grid on;

figure(2)
loglog(odd_var(200),norm_Odd_f_200,'-o',even_var(200),norm_Even_f_200,'-o','MarkerSize',3);
xlabel('m','FontSize',10);
ylabel('N_m', 'FontSize',10);
xticks(10:30:100);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16);
grid on;

figure(3)
loglog(odd_var(199),norm_Odd_f_199,'-o',even_var(199),norm_Even_f_199,'-o','MarkerSize',3);
xlabel('m','FontSize',10);
ylabel('N_m', 'FontSize',10);
xticks(10:30:100);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16);
grid on;

figure(4)
loglog(odd_var(200),norm_Odd_f_dx_200,'-o',even_var(200),norm_Even_f_dx_200,'-o','MarkerSize',3);
xlabel('m','FontSize',10);
ylabel('N_m^{(1)}', 'FontSize',10);
xticks(10:30:100);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16);
grid on;

%% compute the different in the norms
diff_norm_Odd = abs(norm_Odd_f_199 - norm_Odd_f_200(1:length(odd_var(199))))./ ...
                norm_Odd_f_200(1:length(odd_var(199)));
diff_norm_Even = abs(norm_Even_f_199 - norm_Even_f_200(1:length(even_var(199))))./ ...
                norm_Even_f_200(1:length(even_var(199)));

% convert to percent
diff_norm_Odd = diff_norm_Odd * 100;
diff_norm_Even = diff_norm_Even * 100;

figure(5)
loglog(odd_var(199),diff_norm_Odd,'-o',even_var(199),diff_norm_Even,'-o',even_var(199),3,'-v','MarkerSize',3);
%loglog(odd_var(199),norm_Odd_f_199,'-o',odd_var(199),norm_Odd_f_200(1:length(odd_var(199))),'-o','MarkerSize',3);


function f = convg_rate(x_values,y_values)
f = log(y_values(end)/y_values(1))/log(x_values(end)/x_values(1));
end

function id_odd = odd_var(M)
id_odd = 2:2:M;
end

function id_even = even_var(M)
id_even = 1:2:M;
end

% returns the norm of the odd moments and the even moments
function [norm_Odd,norm_Even] = norm_OddEven(data,M,P)
odd_Moments = odd_var(M);
even_Moments = even_var(M);


Odd_Moments = data(odd_Moments+1,:)';
Even_Moments = data(even_Moments+1,:)';

norm_Odd = sqrt(dot(Odd_Moments,P * Odd_Moments,1));
norm_Even = sqrt(dot(Even_Moments,P * Even_Moments,1));
end

function [norm_Odd,norm_Even] = norm_dx_OddEven(data,M,P,D)
odd_Moments = 2 : 2 : M;
even_Moments = 1:2:M;


Odd_Moments = data(odd_Moments+1,:)';
Even_Moments = data(even_Moments+1,:)';

norm_Odd = sqrt(dot(D * Odd_Moments,P * D * Odd_Moments,1));
norm_Even = sqrt(dot(D * Even_Moments,P * D* Even_Moments,1));
end