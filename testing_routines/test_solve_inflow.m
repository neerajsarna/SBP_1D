function test_solve_inflow

par = struct(...
'name','1D Advection',... % name of example
'ic',@ic,... % initial conditions
'relax',@relax,... % production term defined below
'bc_inhomo',@bc_inhomo,... % source term (defined below)
'ax',[0 1],... % coordinates of computational domain
 't_end',0.3,... % the end time of the computation
 'CFL',4.0,...      % the crude cfl number
 'num_bc',2,... % number of boundaries in the domain
 'pres_ID1',true,... % whether we need to prescribe something at x = x_start
 'pres_ID2',true,... % whether we need to prescribe something at x = x_end
 'var_output',1,... % the variable which should be plotted
'output',@output,... % problem-specific output routine (defined below)
'save_during',false, ... % should we save during the computation
'compute_during', @compute_during ...
);

par.t_plot = false;

par.n = 150;

% we need the boundary matrix and the penalty matrix for both the
% boundaries
par.B = cell(par.num_bc,1);
par.penalty = cell(par.num_bc,1);

for k = 3:100
    par.n_eqn = k;  
   
    
    % first we develop the filenames
    filenames = dvlp_filenames(par.n_eqn);
    
    
    % we develop the system data
    [par.Ax,par.B{2}] = get_system_data(filenames);
    
    % stabilise the boundary conditions, was slow in mathematica
    par.B{2} = stabilize_boundary(par.Ax,par.B{2});
    
    % develop the boundary matrix at ID2 or x=0
    par.B{1} = dvlp_B_ID1(par.B{2});
    
    % develop the characteristic penalty matrix
    [par.penalty{1}] = dvlp_penalty(-par.Ax,par.B{1});
    [par.penalty{2}] = dvlp_penalty(par.Ax,par.B{2});
    
    % develop the penalty*B matrix
    for i = 1 : par.num_bc
        par.penalty_B{i} = par.penalty{i} * par.B{i};
    end
    
    result1 = solver(par);
    
    
    par.Ax = dvlp_Ax1D(par.n_eqn);
    par.B{2} = dvlp_BInflow1D(par.n_eqn);
    
    % stabilise the boundary conditions, was slow in mathematica
    par.B{2} = stabilize_boundary(par.Ax,par.B{2});
    
    % develop the boundary matrix at ID2 or x=0
    par.B{1} = dvlp_B_ID1(par.B{2});
    
    % develop the characteristic penalty matrix
    [par.penalty{1}] = dvlp_penalty(-par.Ax,par.B{1});
    [par.penalty{2}] = dvlp_penalty(par.Ax,par.B{2});
    
    % develop the penalty*B matrix
    for i = 1 : par.num_bc
        par.penalty_B{i} = par.penalty{i} * par.B{i};
    end
    
    result2 = solver(par);
    
    diff_result = 0;
    
    for j = 1 : length(result1)
        diff_result = diff_result + norm(result1(j).sol-result2(j).sol);
    end
    
    dlmwrite('error',diff_result,'-append','delimiter','\n');
end

end

% working on inflow boundaries, we consider vacum boundary conditions
function f = bc_inhomo(B,bc_id)
    switch bc_id
        % boundary at x = x_start
        case 1
            thetaIn = 0;
            f = thetaIn * B(:,3)/sqrt(2); 
            
        case 2
            thetaIn = 0;
            f = thetaIn * B(:,3)/sqrt(2); 
    end

end
        
function f = ic(x,id)
if id == 1
    f = exp(-(x-0.5).*(x-0.5)*100);
    
else
    f = zeros(length(x),1);
end
end

function f = relax(x,id)
Kn = 0.1;
f = zeros(length(x),1);

% anything above temperature has to be relaxed
if id > 3
    f = -1/Kn * ones(length(x),1);
end

end

function[filenames] = dvlp_filenames(nEqn)

filenames = struct;

filenames.B = strcat("system_matrices1D/Binflow_1D_",num2str(nEqn));
filenames.B = strcat(filenames.B,".txt");

filenames.Ax = strcat("system_matrices1D/A1_1D_",num2str(nEqn));
filenames.Ax = strcat(filenames.Ax,".txt");
end

% storage of data during time steps becomes expensive, therefore, we only
% focus storing the norms at different times
function compute_during(U,weight,k_RK,PX,DX,t)
norm_f = sqrt(cell2mat(cellfun(@(a) a'*PX*a,U,'Un',0)));
norm_dx_f = sqrt(cell2mat(cellfun(@(a) (DX*a)'*PX*DX*a,U,'Un',0)));

n_eqn = length(U);

norm_dt_f = zeros(1,n_eqn);
% number of RK steps
num_step = length(weight);

% loop over all the components
for i = 1 : n_eqn
    temp = zeros(size(DX,1),1);
    % time derivative of the i-th component at all the points
for j = 1 : num_step
    temp = temp + weight(j) * k_RK{j}{i};
end
    norm_dt_f(i) = sqrt(temp'*PX*temp);
end

filename = strcat('result_Inflow/result_Reference_temp/inflow_t_',num2str(t), ...
                '_points_',num2str(size(DX,2)),'_neqn_',num2str(length(U)),'.txt');

dlmwrite(filename,norm_f,'delimiter','\t','precision',10);
dlmwrite(filename,norm_dx_f,'delimiter','\t','-append','precision',10);
dlmwrite(filename,norm_dt_f,'delimiter','\t','-append','precision',10);
end
