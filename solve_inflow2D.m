function solve_inflow2D(M)

par = struct(...
'name','Inflow2D',... % name of example
'ic',@ic,... % initial conditions
'bc_inhomo',@bc_inhomo,... % source term (defined below)
'ax',[0 1],... % coordinates of computational domain
 't_end',0.3,... % the end time of the computation
 'CFL',2.0,...      % the crude cfl number
 'num_bc',2,... % number of boundaries in the domain
 'pres_ID1',true,... % whether we need to prescribe something at x = x_start
 'pres_ID2',true,... % whether we need to prescribe something at x = x_end
 'var_output',1,... % the variable which should be plotted
'output',@output,... % problem-specific output routine (defined below)
'save_during',false, ... % should we save during the computation
'compute_during', @compute_during, ...
'save_norms', @save_norms, ... % routine for computing and saving the norm
'prod_explicit',true ... % whether we provide an explicit expression for production term or not
);

par.M = M;

if (M < 3)
    error('M should be greater than 3.');
end

% to plot during computation or not to plot, thats the question
par.t_plot = true;

par.n = 50;

if (M == 55)
    par.save_during = true;
end
   
% we need the boundary matrix and the penalty matrix for both the
% boundaries
par.B = cell(par.num_bc,1);
par.penalty = cell(par.num_bc,1);

% first we develop the matrices
par.B{2} = dvlp_BInflow2D(M);
par.Ax = dvlp_Ax2D(M);
par.P = dvlp_Prod2D(M);
par.Kn = inf;

% stabilise the boundary conditions with Onsager
par.B{2} = stabilize_boundary(par.Ax,par.B{2},M);

% develop the boundary matrix at ID2 or x=0
par.B{1} = dvlp_B_ID1(par.B{2},M);

% develop the odd penalty matrix
[par.penalty{2}] = dvlp_penalty_odd(par.Ax,M);
par.penalty{1} = par.penalty{2};

%we remove the variables which are not coupled to density
[par.id_all,par.id_local_odd,par.idx_trun,par.idx_trun_odd,par.idx_trun_even] =  Trun_id(M);
par.Ax = par.Ax(par.id_all,par.id_all);
par.P = par.P(par.id_all,par.id_all);
par.B{1} = par.B{1}(par.id_local_odd,par.id_all);
par.B{2} = par.B{2}(par.id_local_odd,par.id_all);
par.penalty{1} = par.penalty{1}(par.id_all,par.id_local_odd);
par.penalty{2} = par.penalty{2}(par.id_all,par.id_local_odd);

par.n_eqn = size(par.Ax,1);

% develop the penalty*B matrix
for i = 1 : par.num_bc
    par.penalty_B{i} = par.penalty{i} * par.B{i};
end

result = solver(par);

% output_filename = strcat('result_Inflow2D/inflow_tend_', ...
%                         num2str(par.t_end),'_points_',num2str(par.n),'_neqn_');
% output_filename = strcat(output_filename,num2str(M),'.txt');

output_filename = 'result_Inflow2D_Mom.txt';

write_result(result,output_filename);
end

function f = ic(x,id)

if id == 1
    f = exp(-(x-0.5).*(x-0.5)*100);
    %f = zeros(length(x),1);
    
else
    f = zeros(length(x),1);
end

end


%penalty matrix based upon the flux matrix
function [penalty] = dvlp_penalty_odd(Ax,M)

[odd_ID,~] = get_id_Odd(M);

odd_ID = flatten_cell(odd_ID);

% we pick the columns which get multiplied by the odd variables
penalty = Ax(:,odd_ID);
end


% vaccum boundary conditions
function f = bc_inhomo(B,bc_id)
    switch bc_id
        % boundary at x = x_start
        case 1
            thetaIn = 0;
            f = thetaIn * (B(:,3)+B(:,4)+B(:,5))/sqrt(2); 
            
        case 2
            thetaIn = 0;
            f = thetaIn * (B(:,3)+B(:,4)+B(:,5))/sqrt(2); 
    end

end


function [int_f,int_dx_f,int_dt_f] = compute_during(U,weight,k_RK,PX,DX,t,t_Old, ...
                                      idx_trun,idx_trun_odd,idx_trun_even)

% number of elements in idx_trun are M+1                                  
M = length(idx_trun)-1; 
num_step = length(weight);

norm_f = zeros(1,2 * (M+1));
norm_dx_f = zeros(1,2 * (M+1));
norm_dt_f = zeros(1,2 * (M+1));
shift = 0;

for i = 0 : M
    
    if i > 0
        shift = shift + size(idx_trun{i},1);
    end
    
    id_odd = 2*i+2;
    id_even = 2*i + 1;
    
    % for every tensor degree, we first store the norm of the even one and
    % then the odd ones
    
    norm_f(id_even) = sqrt(sum(cell2mat(cellfun(@(a) a'*PX*a,U(shift + idx_trun_even{i+1}),'Un',0))));
    norm_dx_f(id_even) = sqrt(sum(cell2mat(cellfun(@(a) (DX*a)'*PX*DX*a,U(shift + idx_trun_even{i+1}),'Un',0))));
    
    % we now loop over all the even variables
    for j = 1 : length(idx_trun_even{i+1})
        
        temp = zeros(size(DX,1),1);
        
        for k = 1 : num_step
            temp = temp + weight(k) * k_RK{k}{shift + idx_trun_even{i+1}(j)};
        end
        
        norm_dt_f(id_even) = norm_dt_f(id_even) + temp'*PX*temp;
    end
    
    norm_dt_f(id_even) = sqrt(norm_dt_f(id_even));
    
    % we first loop over all the odd variables
    if i > 0
        norm_dx_f(id_odd) = sqrt(sum(cell2mat(cellfun(@(a) (DX*a)'*PX*DX*a,U(shift + idx_trun_odd{i+1}),'Un',0))));
        norm_f(id_odd) = sqrt(sum(cell2mat(cellfun(@(a) a'*PX*a,U(shift + idx_trun_odd{i+1}),'Un',0))));
        
        for j = 1 : length(idx_trun_odd{i+1})
            temp = zeros(size(DX,1),1);
            
            for k = 1 : num_step
                temp = temp + weight(k) * k_RK{k}{shift + idx_trun_odd{i+1}(j)};
            end
            
            norm_dt_f(id_odd) = norm_dt_f(id_odd) + temp'*PX*temp;
        end
    end
    
    norm_dt_f(id_odd) = sqrt(norm_dt_f(id_odd));
   
    
end

% we perform the integral here itself 
% using the midpoint rule
int_f = norm_f * (t-t_Old);
int_dx_f = norm_dx_f * (t-t_Old);
int_dt_f = norm_dt_f * (t-t_Old);

end

function save_norms(int_f,int_dx_f,int_dt_f,n,M)
filename = strcat('result_Inflow2D/result_Reference/inflow_norms', ...
                '_points_',num2str(n),'_neqn_',num2str(M),'.txt');

dlmwrite(filename,int_f,'delimiter','\t','precision',10);
dlmwrite(filename,int_dx_f,'delimiter','\t','-append','precision',10);
dlmwrite(filename,int_dt_f,'delimiter','\t','-append','precision',10);
end

% for the given initial and boundary conditions, we can truncate the
% matrices to make them smaller
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
