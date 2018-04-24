% we compute the expected rate of convergence
% details of the input argumens 
% 1. n_eqn = reference number of equations
% 2. t_end = end time at which we do the computations
% 3. n = number of grid points along the spatial domain
% 4. foldername = folder in which the results are
function [expected_rate,loc_truncate] = expected_convg_rate(n_eqn,n_eqn2,t_end,n,foldername,filename)


% filename of the file which contains all the norms 
norms_filename = strcat(foldername,'/result_Reference/',filename,'_norms_points_',num2str(n),'_neqn_',...
                            num2str(n_eqn),'.txt');

result = dlmread(norms_filename,'\t');

int_f = result(1,:);
int_dx_f = result(2,:);
int_dt_f = result(3,:);

delta_x = 1/n;

% % we need the P matrix for computing the error
[D,P,~] = sbp_traditional_2(delta_x,n);


%% we also extract the reference solution at t = 0.3
output_filename = strcat(foldername,'/',filename,'_tend_',num2str(t_end),'_points_',num2str(n),'_neqn_');
output_filename = strcat(output_filename,num2str(n_eqn),'.txt');
    
result_final = dlmread(output_filename,'\t');
result_final = result_final(2:end,:)';

% solution corresponding to the second reference
output_filename = strcat(foldername,'/',filename,'_tend_',num2str(t_end),'_points_',num2str(n),'_neqn_');
output_filename = strcat(output_filename,num2str(n_eqn2),'.txt');

result_final_2 = dlmread(output_filename,'\t');
result_final_2 = result_final_2(2:end,:)';

norm_final = sqrt(dot(result_final,P*result_final,1));

% 
% % plot the variation in norm f
% 
loc_truncate = find_truncate(result_final,result_final_2,P,D);

id_odd = odd_var(loc_truncate);
id_even = even_var(loc_truncate);

id_odd_full = odd_var(n_eqn);
id_even_full = even_var(n_eqn);

%% polyfit for f
[P_Odd_f,y_Odd_f] = polyfit_linear(log(id_odd),log(int_f(id_odd)));
[P_Even_f,y_Even_f] = polyfit_linear(log(id_even),log(int_f(id_even)));
% 

%% polyfit for dxf
[P_Odd_dx_f,y_Odd_dx_f] = polyfit_linear(log(id_odd),log(int_dx_f(id_odd)));
[P_Even_dx_f,y_Even_dx_f] = polyfit_linear(log(id_even),log(int_dx_f(id_even)));

% 
%% polyfit for dt_f
[P_Odd_dt_f,y_Odd_dt_f] = polyfit_linear(log(id_odd),log(int_dt_f(id_odd)));
[P_Even_dt_f,y_Even_dt_f] = polyfit_linear(log(id_even),log(int_dt_f(id_even)));

%% polyfit for f at t = t_end
[P_Odd_f_final,y_Odd_f_final] = polyfit_linear(log(id_odd),log(norm_final(id_odd)));
[P_Even_f_final,y_Even_f_final] = polyfit_linear(log(id_even),log(norm_final(id_even)));


% 
%% plotting for f
figure(1)
loglog(id_odd_full,int_f(id_odd_full),'-o',id_even_full,int_f(id_even_full),'-o',id_odd,exp(y_Odd_f),'r-*',...
                id_even,exp(y_Even_f),'k-*',...
                    'MarkerSize',3);
legend('odd moments','even moments','linear fit odd','linear fit even');
xlim([1 n_eqn]);
title('norm int f');
grid on;
% 
%% plotting for dt f
figure(2)
loglog(id_odd_full,int_dt_f(id_odd_full),'-o',id_even_full,int_dt_f(id_even_full),'-o',id_odd,exp(y_Odd_dt_f),'r-*',...
                id_even,exp(y_Even_dt_f),'k-*',...
                    'MarkerSize',3);
legend('odd moments','even moments','linear fit odd','linear fit even');
xlim([1 n_eqn]);
title('norm int dt f');
grid on;
% 
 %% plotting for dx f
figure(3)
loglog(id_odd_full,int_dx_f(id_odd_full),'-o',id_even_full,int_dx_f(id_even_full),'-o',id_odd,exp(y_Odd_dx_f),'r-*',...
                id_even,exp(y_Even_dx_f),'k-*',...
                    'MarkerSize',3);
legend('odd moments','even moments','linear fit odd','linear fit even');
xlim([1 n_eqn]);
title('norm int dx f');
grid on;
% 

%% plotting for f at t = t_end
figure(4)
loglog(id_odd_full,norm_final(id_odd_full),'-o',id_even_full,norm_final(id_even_full),'-o',id_odd,exp(y_Odd_f_final),'r-*',...
                id_even,exp(y_Even_f_final),'k-*',...
                    'MarkerSize',3);
legend('odd moments','even moments','linear fit odd','linear fit even');
xlim([1 n_eqn]);
title('norm f(t = t_{end})');
grid on;


%% display the convergence rates 
quantities = {'int(f)';'int(d_xf)';'int(d_tf)';'f(t_{end})'};
even_Order = abs([P_Even_f(1), P_Even_dx_f(1),P_Even_dt_f(1),P_Even_f_final(1)])';
odd_Order = abs([P_Odd_f(1), P_Odd_dx_f(1),P_Odd_dt_f(1),P_Odd_f_final(1)])';
reduction = even_Order - [0.5,1,0.5,0.5]';
expected_rate = min(reduction);

T = table(quantities,even_Order,odd_Order,reduction);
disp(T);
disp('EXPECTED RATE');
disp(expected_rate);

end

function id_odd = odd_var(M)
id_odd = 2:2:M;
end

function id_even = even_var(M)
id_even = 1:2:M;
end

% returns the polyfit with a linear fit and the corresponding y function
function [P,yfit] = polyfit_linear(x,y)
P = polyfit(x,y,1);
yfit = P(1) * x + P(2);
end

function f = find_truncate(result,result_2,P,D)

cutoff_per = 5;

norm1 = sqrt(dot(result,P*result,1));
norm2 = sqrt(dot(result_2,P*result_2,1));

norm1_dx = sqrt(dot(D * result,P* D * result,1));
norm2_dx = sqrt(dot(D * result_2,P* D * result_2,1));

min_dim = min(length(norm1),length(norm2));
diff1 = abs(norm1(1:min_dim)-norm2(1:min_dim))*100./norm1(1:min_dim);
diff1_dx = abs(norm1_dx(1:min_dim)-norm2_dx(1:min_dim))*100./norm1_dx(1:min_dim);

f = min(find(diff1 > cutoff_per,1),find(diff1_dx > cutoff_per,1));
end

