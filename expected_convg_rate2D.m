% we compute the expected rate of convergence
% details of the input argumens 
% 1. n_eqn = reference number of equations
% 2. t_end = end time at which we do the computations
% 3. n = number of grid points along the spatial domain
% 4. foldername = folder in which the results are
function [expected_rate,loc_truncate] = expected_convg_rate2D(M,t_end,n,foldername,filename)


% filename of the file which contains all the norms 
norms_filename = strcat(foldername,'/result_Reference/',filename,'_norms_points_',num2str(n),'_neqn_',...
                            num2str(M),'.txt');

result = dlmread(norms_filename,'\t');

int_f = result(1,:);
int_dx_f = result(2,:);
int_dt_f = result(3,:);

delta_x = 1/n;

% % we need the P matrix for computing the error
[~,P,~] = sbp_traditional_2(delta_x,n);


%% we also extract the reference solution at t = 0.3
output_filename = strcat(foldername,'/',filename,'_tend_',num2str(t_end),'_points_',num2str(n),'_neqn_');
output_filename = strcat(output_filename,num2str(M),'.txt');
    
result_final = dlmread(output_filename,'\t');
result_final = result_final(2:end,:)';


[~,~,idx_trun,idx_trun_odd,idx_trun_even] =  Trun_id(M);

norm_final = zeros(2*M+2,1);
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
    else
            norm_final(id_even) = sqrt(sum(dot(result_final(:,shift + idx_trun_even{i+1}), ...
                        P * result_final(:,shift + idx_trun_even{i+1}),1)));
    end
    
    
    % we first loop over all the odd variables
    if isempty(idx_trun_odd{i+1})
      norm_final(id_odd) = 0;
    else
     norm_final(id_odd) = sqrt(sum(dot(result_final(:,shift + idx_trun_odd{i+1}), ...
                               P * result_final(:,shift + idx_trun_odd{i+1}),1)));
    end
    
end

 
id_odd_full = odd_var(2 * M + 2);
id_even_full = even_var(2 * M + 2);

id_odd = id_odd_full;
id_even = id_even_full;

%% polyfit for f
[P_Odd_f,y_Odd_f] = polyfit_linear(log(id_odd),log(int_f(id_odd)));
[P_Even_f,y_Even_f] = polyfit_linear(log(id_even),log(int_f(id_even)));
% 

%% polyfit for dxf
[P_Odd_dx_f,y_Odd_dx_f] = polyfit_linear(log(id_odd),log(int_dx_f(id_odd)));
[P_Even_dx_f,y_Even_dx_f] = polyfit_linear(log(id_even),log(int_dx_f(id_even)));

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
xlim([0 M+2]);
title('norm int f');
grid on;
% 
%% plotting for dt f
figure(2)
loglog(id_odd_full,int_dt_f(id_odd_full),'-o',id_even_full,int_dt_f(id_even_full),'-o',id_odd,exp(y_Odd_dt_f),'r-*',...
                id_even,exp(y_Even_dt_f),'k-*',...
                    'MarkerSize',3);
legend('odd moments','even moments','linear fit odd','linear fit even');
xlim([0 M+2]);
title('norm int dt f');
grid on;
% 
 %% plotting for dx f
figure(3)
loglog(id_odd_full,int_dx_f(id_odd_full),'-o',id_even_full,int_dx_f(id_even_full),'-o',id_odd,exp(y_Odd_dx_f),'r-*',...
                id_even,exp(y_Even_dx_f),'k-*',...
                    'MarkerSize',3);
legend('odd moments','even moments','linear fit odd','linear fit even');
xlim([0 M+2]);
title('norm int dx f');
grid on;
% 

%% plotting for f at t = t_end
figure(4)
loglog(id_odd_full,norm_final(id_odd_full),'-o',id_even_full,norm_final(id_even_full),'-o',id_odd,exp(y_Odd_f_final),'r-*',...
                id_even,exp(y_Even_f_final),'k-*',...
                    'MarkerSize',3);
legend('odd moments','even moments','linear fit odd','linear fit even');
xlim([0 M+2]);
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

