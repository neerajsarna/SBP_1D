% we compute the expected rate of convergence
% details of the input argumens 
% 1. n_eqn = reference number of equations
% 2. t_end = end time at which we do the computations
% 3. n = number of grid points along the spatial domain
% 4. foldername = folder in which the results are
function [expected_rate,loc_truncate] = expected_convg_rate2D(M,t_end,n,...
                                            foldername,filename)


% 


% filename of the file which contains all the norms 
norms_filename = strcat(foldername,'/result_Reference/',filename,'_norms_points_',num2str(n),'_neqn_',...
                            num2str(M),'.txt');

result = dlmread(norms_filename,'\t');

int_f = result(1,:);
int_dx_f = result(2,:);
int_dt_f = result(3,:);

delta_x = 1/n;

% % we need the P matrix for computing the error
[DX,P,~] = sbp_traditional_2(delta_x,n);


%% we also extract the reference solution at t = 0.3
output_filename = strcat(foldername,'/',filename,'_tend_',num2str(t_end),'_points_',num2str(n),'_neqn_');
output_filename = strcat(output_filename,num2str(M),'.txt');
    
result_final = dlmread(output_filename,'\t');
result_final = result_final(2:end,:)';


[~,~,idx_trun,idx_trun_odd,idx_trun_even] =  Trun_id(M);

norm_final = zeros(1,2*M+2);
norm_dx_final = zeros(1,2*M+2);

shift = 0;
for i = 0 : M
    
    if i > 0
        shift = shift + size(idx_trun{i},1);
    end
    
    id_odd = 2*i+2;
    id_even = 2*i + 1;
    
    % for every tensor degree, we first store the norm of the even one and
    % then the odd ones
    
    if isempty(idx_trun_even{i+1})
            norm_final(id_even) = 0;
            norm_dx_final(id_even) = 0;
    else
        norm_final(id_even) = sqrt(sum(dot(result_final(:,shift + idx_trun_even{i+1}), ...
            P * result_final(:,shift + idx_trun_even{i+1}),1)));
        
        norm_dx_final(id_even) = sqrt(sum(dot(DX * result_final(:,shift + idx_trun_even{i+1}), ...
            P * DX * result_final(:,shift + idx_trun_even{i+1}),1)));
    end
    
    
    % we first loop over all the odd variables
    if isempty(idx_trun_odd{i+1})
      norm_final(id_odd) = 0;
      norm_dx_final(id_odd) = 0;
    else
     norm_final(id_odd) = sqrt(sum(dot(result_final(:,shift + idx_trun_odd{i+1}), ...
                               P * result_final(:,shift + idx_trun_odd{i+1}),1)));
                           
     norm_dx_final(id_odd) = sqrt(sum(dot(DX * result_final(:,shift + idx_trun_odd{i+1}), ...
                               P * DX * result_final(:,shift + idx_trun_odd{i+1}),1)));                      
    end
    
end

% given a tensor degree, we either have odd or even moment. 
% find the locations where non zero values of magnitude 
id_odd = odd_var(2 * M + 2);
id_even = even_var(2 * M + 2);

%% extract the norms which are not zero
loc_truncate = M+1;

[M_values_Odd,int_f_Odd] = return_nz_Odd(id_odd,int_f(id_odd),loc_truncate);
[M_values_Even,int_f_Even] = return_nz_Even(id_even,int_f(id_even),loc_truncate);

[~,int_dx_f_Odd] = return_nz_Odd(id_odd,int_dx_f(id_odd),loc_truncate);
[~,int_dx_f_Even] = return_nz_Even(id_even,int_dx_f(id_even),loc_truncate);

[~,int_dt_f_Odd] = return_nz_Odd(id_odd,int_dt_f(id_odd),loc_truncate);
[~,int_dt_f_Even] = return_nz_Even(id_even,int_dt_f(id_even),loc_truncate);

[~,f_Odd] = return_nz_Odd(id_odd,norm_final(id_odd),loc_truncate);
[~,f_Even] = return_nz_Even(id_even,norm_final(id_even),loc_truncate);

[~,dx_f_Odd] = return_nz_Odd(id_odd,norm_dx_final(id_odd),loc_truncate);
[~,dx_f_Even] = return_nz_Even(id_even,norm_dx_final(id_even),loc_truncate);

%% find the location where we would like to truncate
loc_truncate = M-10;
loc_truncate_Odd = find(M_values_Odd >= loc_truncate,1);
loc_truncate_Even = find(M_values_Even >= loc_truncate,1);

%% polyfit for f,
[P_Odd_f,y_Odd_f] = polyfit_linear(log(M_values_Odd(4:loc_truncate_Odd)), ...
                                    log(int_f_Odd(4:loc_truncate_Odd)));
                                
[P_Even_f,y_Even_f] = polyfit_linear(log(M_values_Even(4:loc_truncate_Even)), ...
                                    log(int_f_Even(4:loc_truncate_Even)));

%% polyfit for dxf
[P_Odd_dx_f,y_Odd_dx_f] = polyfit_linear(log(M_values_Odd(3:loc_truncate_Odd)), ...
                                        log(int_dx_f_Odd(3:loc_truncate_Odd)));
                                    
[P_Even_dx_f,y_Even_dx_f] = polyfit_linear(log(M_values_Even(3:loc_truncate_Even)), ...
                                        log(int_dx_f_Even(3:loc_truncate_Even)));

%% polyfit for dt_f
[P_Odd_dt_f,y_Odd_dt_f] = polyfit_linear(log(M_values_Odd(3:loc_truncate_Odd)), ...
                                        log(int_dt_f_Odd(3:loc_truncate_Odd)));
                                    
[P_Even_dt_f,y_Even_dt_f] = polyfit_linear(log(M_values_Even(3:loc_truncate_Even)), ...
                                            log(int_dt_f_Even(3:loc_truncate_Even)));

%% polyfit for f at t = t_end
[P_Odd_f_final,y_Odd_f_final] = polyfit_linear(log(M_values_Odd(3:loc_truncate_Odd)), ...
                                                log(f_Odd(3:loc_truncate_Odd)));
                                            
[P_Even_f_final,y_Even_f_final] = polyfit_linear(log(M_values_Even(4:loc_truncate_Even)), ...
                                                log(f_Even(4:loc_truncate_Even)));

% polyfit for dx_f at t = t_end
[P_Odd_dx_f_final,y_Odd_dx_f_final] = polyfit_linear(log(M_values_Odd(3:loc_truncate_Odd)), ...
                                                log(dx_f_Odd(3:loc_truncate_Odd)));
                                            
[P_Even_dx_f_final,y_Even_dx_f_final] = polyfit_linear(log(M_values_Even(4:loc_truncate_Even)), ...
                                                log(dx_f_Even(4:loc_truncate_Even)));
% 
%% plotting for f
figure(1)
loglog(M_values_Odd,int_f_Odd,'-o',M_values_Even,int_f_Even,'-o', ...
                M_values_Odd(4:loc_truncate_Odd),exp(y_Odd_f),'r-*',...
                M_values_Even(4:loc_truncate_Even),exp(y_Even_f),'k-*',...
                    'MarkerSize',3);
                
legend('odd moments','even moments','linear fit odd', ...
        'linear fit even','Location','best');
xlim([0 M+2]);
title('variation of N_m');
xlabel('m+1');
ylabel('N_m');
xt = get(gca, 'YTick');
set(gca, 'FontSize', 16);
grid on;

% 
%% plotting for dx f
figure(2)
loglog(M_values_Odd,int_dx_f_Odd,'-o',M_values_Even,int_dx_f_Even,'-o', ...
            M_values_Odd(3:loc_truncate_Odd),exp(y_Odd_dx_f),'r-*',...
                M_values_Even(3:loc_truncate_Even),exp(y_Even_dx_f),'k-*',...
                    'MarkerSize',3);
legend('odd moments','even moments','linear fit odd', ...
        'linear fit even','Location','best');
xlim([0 M+2]);
title('variation of N_m^{(x)}');
xlabel('m+1');
ylabel('N_m^{(x)}');
xt = get(gca, 'YTick');
set(gca, 'FontSize', 16);
grid on;
% 
 %% plotting for dt f
figure(3)
loglog(M_values_Odd,int_dt_f_Odd,'-o',M_values_Even,int_dt_f_Even,'-o', ...
            M_values_Odd(3:loc_truncate_Odd),exp(y_Odd_dt_f),'r-*',...
                M_values_Even(3:loc_truncate_Even),exp(y_Even_dt_f),'k-*',...
                    'MarkerSize',3);
legend('odd moments','even moments','linear fit odd', ...
        'linear fit even','Location','best');
xlim([0 M+2]);
title('variation of N_m^{(t)}');
xlabel('m+1');
ylabel('N_m^{(t)}');
xt = get(gca, 'YTick');
set(gca, 'FontSize', 16);
grid on;
% 

%% plotting for f at t = t_end
figure(4)
loglog(M_values_Odd,f_Odd,'-o',M_values_Even,f_Even,'-o', ...
                M_values_Odd(3:loc_truncate_Odd),exp(y_Odd_f_final),'r-*',...
                M_values_Even(4:loc_truncate_Even),exp(y_Even_f_final),'k-*',...
                'MarkerSize',3);
            
legend('odd moments','even moments','linear fit odd', ...
        'linear fit even','Location','best');
xlim([0 M+2]);
title('variation of N_m^{(T)}');
xlabel('m+1');
ylabel('N_m^{(T)}');
xt = get(gca, 'YTick');
set(gca, 'FontSize', 16);
grid on;


%% display the convergence rates 
quantities = {'int(f)';'int(d_xf)';'int(d_tf)';'f(t_{end})'};
even_Order = abs([P_Even_f(1), P_Even_dx_f(1),P_Even_dt_f(1),P_Even_f_final(1)])';
odd_Order = abs([P_Odd_f(1), P_Odd_dx_f(1),P_Odd_dt_f(1),P_Odd_f_final(1)])';
reduction_Even = even_Order - [0.5,1,0.5,0.5]';
reduction_Odd = odd_Order - [0.5,0.5,0.5,0.5]';
expected_rate = min([reduction_Even;reduction_Odd(2)]);

T = table(quantities,even_Order,odd_Order,reduction_Even,reduction_Odd);
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


function [id_all,id_local_odd,idx_trun,idx_trun_local_odd,idx_trun_local_even] =  Trun_id(M)
all_idx = cell(M+1,1);
idx_odd = cell(M+1,1);

idx_trun = cell(M+1,1);
idx_trun_local_odd = cell(M+1,1);
idx_trun_local_even = cell(M+1,1);

%% we collect the indices of the basis functions
for i = 0 : M
    [all_idx{i+1},~,~] = IDX_Full(i);
    idx_odd{i+1} = all_idx{i+1}(rem(all_idx{i+1}(:,1),2) ~= 0,:);
end

id_all = [];
id_local_odd = [];

%% we now truncate
for i = 0 : M
    shift_id = sum(cell2mat(cellfun(@(a) size(a,1),all_idx(1:i),'Un',0)));
    
    % all the basis which have 0 or 2 or the permutation of the two in y
    % and z direction will contribute
    [~,ia,~] = intersect(all_idx{i+1}(:,2:end),[0, 0;2,0;0,2],'rows');
    
    % location of the variables
    idx_trun{i + 1} = all_idx{i+1}(sort(ia),:);
    idx_trun_local_odd{i + 1} = find(rem(idx_trun{i+1}(:,1),2) ~= 0);
    idx_trun_local_even{i + 1} = find(rem(idx_trun{i+1}(:,1),2) == 0);
    
    
    id_all = [id_all sort(ia)'+shift_id];
    
    % same for the odd variables
    shift_id = sum(cell2mat(cellfun(@(a) size(a,1),idx_odd(1:i),'Un',0)));
    
    [~,ia,~] = intersect(idx_odd{i+1}(:,2:end),[0, 0;2,0;0,2],'rows');
    
    id_local_odd = [id_local_odd sort(ia)'+shift_id];
end 

end

% we return the M_values of non zero odd norms
function [M_values,norm_values] = return_nz_Odd(id_Odd,norm_Odd,loc_truncate)

loc_nz = find(norm_Odd ~= 0);
id_Odd_nz = id_Odd(loc_nz);

M_values = id_Odd_nz/2; 
norm_values = norm_Odd(loc_nz);

% only keep values which are less than loc_truncate
loc_less_truncate = find(M_values <= loc_truncate);

M_values = M_values(loc_less_truncate);
norm_values = norm_values(loc_less_truncate);

end

% return the M values of non zero even norms
function [M_values,norm_values] = return_nz_Even(id_Even,norm_Even,loc_truncate)

loc_nz = find(norm_Even ~= 0);
id_Even_nz = id_Even(loc_nz);

M_values = ceil(id_Even_nz/2); 
norm_values = norm_Even(loc_nz);

% only keep values which are less than loc_truncate
loc_less_truncate = find(M_values <= loc_truncate);

M_values = M_values(loc_less_truncate);
norm_values = norm_values(loc_less_truncate);
end
