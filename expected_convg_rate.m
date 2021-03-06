% we compute the expected rate of convergence
% details of the input argumens 
% 1. n_eqn = reference number of equations
% 2. t_end = end time at which we do the computations
% 3. n = number of grid points along the spatial domain
% 4. foldername = folder in which the results are
function [expected_rate] = expected_convg_rate(n_eqn,t_end,n, ...
                                        foldername,filename)


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
output_filename = strcat(foldername,'/',filename,'_tend_', ...
                        num2str(t_end),'_points_',num2str(n),'_neqn_');
output_filename = strcat(output_filename,num2str(n_eqn),'.txt');
    
result_final = dlmread(output_filename,'\t');
result_final = result_final(2:end,:)';

norm_final = sqrt(dot(result_final,P*result_final,1));
             

id_odd = odd_var(n_eqn);
id_even = even_var(n_eqn);

%% the points to be used during the linear fit, we start from three to ignore the initial part
id_full_E = 3:1:180;
id_full_O = 3:1:100;

id_odd_full = find(mod(id_full_O,2) == 0);
id_even_full = find(mod(id_full_E,2) ~=0 );

%% polyfit for f

[P_f_O,y_f_O] = polyfit_linear(log(id_full_O(id_odd_full)), ...
                                log(int_f(id_full_O(id_odd_full))));

[P_f_E,y_f_E] = polyfit_linear(log(id_full_E(id_even_full)), ...
                                log(int_f(id_full_E(id_even_full))));

% 

%% polyfit for dxf
[P_dx_f_O,y_dx_f_O] = polyfit_linear(log(id_full_O(id_odd_full)), ...
                                    log(int_dx_f(id_full_O(id_odd_full))));
[P_dx_f_E,y_dx_f_E] = polyfit_linear(log(id_full_E(id_even_full)), ...
                                    log(int_dx_f(id_full_E(id_even_full))));
% 
%% polyfit for dt_f
[P_dt_f_O,y_dt_f_O] = polyfit_linear(log(id_full_O(id_odd_full)), ...
                                    log(int_dt_f(id_full_O(id_odd_full))));
                                
[P_dt_f_E,y_dt_f_E] = polyfit_linear(log(id_full_E(id_even_full)), ...
                                    log(int_dt_f(id_full_E(id_even_full))));
%% polyfit for f at t = t_end
[P_f_final_O,y_f_final_O] = polyfit_linear(log(id_full_O(id_odd_full)), ...
                                        log(norm_final(id_full_O(id_odd_full))));
[P_f_final_E,y_f_final_E] = polyfit_linear(log(id_full_E(id_even_full)), ...
                                            log(norm_final(id_full_E(id_even_full))));

%% plotting for f
figure(1)
fig = loglog(id_odd,int_f(id_odd),'-o', ...
             id_even,int_f(id_even),'-square',...
             id_full_O(id_odd_full),exp(y_f_O),'-*',...
             id_full_E(id_even_full),exp(y_f_E),'->',...
             'MarkerSize',4);
                
legend('odd moments','even moments','linear fit odd','linear fit even');
xlim([1 n_eqn]);
title('variation of N_m');
xlabel('m+1');
ylabel('N_m');
xt = get(gca, 'YTick');
set(gca, 'FontSize', 16);
grid on;


% 
%% plotting for dt f
figure(2)
fig = loglog(id_odd,int_dt_f(id_odd),'-o', ...
             id_even,int_dt_f(id_even),'-square',...
                id_full_O(id_odd_full),exp(y_dt_f_O),'-*',...
                id_full_E(id_even_full),exp(y_dt_f_E),'->',...
                    'MarkerSize',4);
legend('odd moments','even moments','linear fit odd','linear fit even');
xlim([1 n_eqn]);
title('variation of N_m^{(t)}');
xlabel('m+1');
ylabel('N_m^{(t)}');
xt = get(gca, 'YTick');
set(gca, 'FontSize', 16);
grid on;

% 

 %% plotting for dx f
figure(3)
fig = loglog(id_odd,int_dx_f(id_odd),'-o', ...
             id_even,int_dx_f(id_even),'-square',...
              id_full_O(id_odd_full),exp(y_dx_f_O),'-*',...
              id_full_E(id_even_full),exp(y_dx_f_E),'->',...
               'MarkerSize',4);
legend('odd moments','even moments','linear fit Odd','linear fit Even');
xlim([1 n_eqn]);
title('variation of N_m^{(x)}');
xlabel('m+1');
ylabel('N_m^{(x)}');
xt = get(gca, 'YTick');
set(gca, 'FontSize', 16);
grid on;

% 

%% plotting for f at t = t_end
figure(4)
fig = loglog(id_odd,norm_final(id_odd),'-o', ...
            id_even,norm_final(id_even),'-square',...
              id_full_O(id_odd_full),exp(y_f_final_O),'-*',...
              id_full_E(id_even_full),exp(y_f_final_E),'->',...
                    'MarkerSize',4);
legend('odd moments','even moments','linear fit odd','linear fit even');
xlim([1 n_eqn]);
title('variation of N_m^{(T)}');
xlabel('m+1');
ylabel('N_m^{(T)}');
xt = get(gca, 'YTick');
set(gca, 'FontSize', 16);
grid on;



%% display the convergence rates 
quantities = {'int(f)';'int(d_xf)';'int(d_tf)';'f(t_{end})'};
odd_Order = abs([P_f_O(1), P_dx_f_O(1),P_dt_f_O(1),P_f_final_O(1)])';
even_Order = abs([P_f_E(1), P_dx_f_E(1),P_dt_f_E(1),P_f_final_E(1)])';
reduction_full = [even_Order(1),even_Order(2),min(odd_Order(3),even_Order(3)),min(odd_Order(4),even_Order(4))]'...
                - [0.5,1,0.5,0.5]';
expected_rate = min(reduction_full);

T = table(quantities,odd_Order,even_Order,reduction_full);
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

