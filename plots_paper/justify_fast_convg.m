clear all;

M_mom_values = 3:2:27;

norm_values = [];

M_ref = 55;

result_M24 = dlmread(strcat('result_HC2D/hc_tend_1_points_300_neqn_',num2str(M_ref), ...
                    '.txt'));

[id_all_M24,id_local_odd_M24,idx_global_even_M24,idx_global_odd_M24, ...
    idx_trun_local_M24,idx_trun_local_odd_M24,idx_trun_local_even_M24] =  Trun_id(M_ref);
B_M24 = dvlp_BWall2D(M_ref);
Boe_M24 = -B_M24(id_local_odd_M24,idx_global_even_M24);
% extract the Aoe matrix
Ax = dvlp_Ax2D(M_ref);
Aoe = Ax(idx_global_odd_M24,idx_global_even_M24);
% Through the present file we try to justify the fast convergence seen for
% the heat conduction problem. 
for M_mom = M_mom_values


[id_all_M3,id_local_odd_M3,idx_global_even_M3,~, ...
    idx_trun_local_M3,idx_trun_local_odd_M3,idx_trun_local_even_M3] =  Trun_id(M_mom);

B_M3 = dvlp_BWall2D(M_mom);

Boe_M3 = -B_M3(id_local_odd_M3,idx_global_even_M3);

% part of Boe_M24 not in Boe_M3
Boe_M24_minus_M3 = Boe_M24(1:end,size(Boe_M3,2)+1:end);

% multiply the remaining part of Boe with the remaining of the even moments
even_result = extract_even_result(result_M24,idx_trun_local_M24, ...
                                    idx_trun_local_even_M24,M_ref);

extra_even = size(Boe_M24,2)-size(Boe_M3,2);
Boe_mom = Boe_M24_minus_M3 * ...
            even_result(end-(extra_even-1):end)';

norm_values = [norm_values norm(Boe_mom)];

disp(M_mom);

end

loglog(M_mom_values,norm_values);


function [id_all,id_local_odd,id_global_even,id_global_odd, ...
        idx_trun_local,idx_trun_local_odd,idx_trun_local_even] =  Trun_id(M)

all_idx = cell(M+1,1);
idx_odd = cell(M+1,1);
idx_even = cell(M+1,1);


idx_trun_local = cell(M+1,1);
idx_trun_local_odd = cell(M+1,1);
idx_trun_local_even = cell(M+1,1);

%% we collect the indices of the basis functions
for i = 0 : M
    [all_idx{i+1},~,~] = IDX_Full(i);
    idx_odd{i+1} = all_idx{i+1}(rem(all_idx{i+1}(:,1),2) ~= 0,:);
    idx_even{i+1} = all_idx{i+1}(rem(all_idx{i+1}(:,1),2) == 0,:);
end

id_all = [];
id_local_odd = [];
id_global_even = [];
id_global_odd = [];

%% we now truncate
for i = 0 : M
    % will get not shifted for i = 0
    shift_id = sum(cell2mat(cellfun(@(a) size(a,1),all_idx(1:i),'Un',0)));
    
    % all the basis which have 0 or 2 or the permutation of the two in y
    % and z direction will contribute
    [~,ia,~] = intersect(all_idx{i+1}(:,2:end),[0, 0;2,0;0,2],'rows');
   
    
    id_all = [id_all sort(ia)'+shift_id];
    
    % location of the variables without the shift 
    idx_trun_local{i + 1} = all_idx{i+1}(sort(ia),:);
    idx_trun_local_odd{i + 1} = find(rem(idx_trun_local{i+1}(:,1),2) ~= 0);
    idx_trun_local_even{i + 1} = find(rem(idx_trun_local{i+1}(:,1),2) == 0);
    
    % same for the odd variables
    shift_id = sum(cell2mat(cellfun(@(a) size(a,1),idx_odd(1:i),'Un',0)));
    
    [~,ia,~] = intersect(idx_odd{i+1}(:,2:end),[0, 0;2,0;0,2],'rows');
    
    % location of variables with the shift
    id_local_odd = [id_local_odd sort(ia)'+shift_id];
    
    shift_id = sum(cell2mat(cellfun(@(a) size(a,1),all_idx(1:i),'Un',0)));
    
    [~,ia,~] = intersect(idx_even{i+1}(:,2:end),[0, 0;2,0;0,2],'rows');
    
    id_global_even = [id_global_even sort(ia)'+shift_id];
    
    [~,ia,~] = intersect(idx_odd{i+1}(:,2:end),[0, 0;2,0;0,2],'rows');
    
    id_global_odd = [id_global_odd sort(ia)'+shift_id];
end 

end

% given the result, find the value of the even moments
function f = extract_even_result(result_M,idx_trun_local,idx_trun_local_even,M)

f = [];
shift = 0;

for i = 0 : M
    
    if i > 0
        shift = shift + size(idx_trun_local{i},1);
    end
        
    % for every tensor degree, we first store the norm of the even one and
    % then the odd ones
    
    f = [f result_M(shift+idx_trun_local_even{i+1},1)'];
    
    
end


end

