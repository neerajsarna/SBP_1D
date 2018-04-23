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
 'var_output',2,... % the variable which should be plotted
'output',@output,... % problem-specific output routine (defined below)
'save_during',true, ... % should we save during the computation
'compute_during', @compute_during, ...
'prod_explicit',true ... % whether we provide an explicit expression for production term or not
);

if (M < 3)
    error('M should be greater than 3');
end
% we don't plot during the computation
par.t_plot = false;

par.n = 300;
   
% we need the boundary matrix and the penalty matrix for both the
% boundaries
par.B = cell(par.num_bc,1);
par.penalty = cell(par.num_bc,1);

% first we develop the matrices
par.B{2} = dvlp_BInflow2D(M);
par.Ax = dvlp_Ax2D(M);
par.P = dvlp_Prod2D(M);
par.Kn = 0.1;

par.n_eqn = size(par.Ax,1);

% stabilise the boundary conditions with Onsager
par.B{2} = stabilize_boundary(par.Ax,par.B{2},M);

% develop the boundary matrix at ID2 or x=0
par.B{1} = dvlp_B_ID1(par.B{2},M);

% develop the odd penalty matrix
[par.penalty{1}] = dvlp_penalty_odd(par.Ax,M);
par.penalty{2} = par.penalty{1};

% develop the penalty*B matrix
for i = 1 : par.num_bc
    par.penalty_B{i} = par.penalty{i} * par.B{i};
end

result = solver(par);

output_filename = strcat('result_Inflow2D/inflow_tend_',num2str(par.t_end),'_points_',num2str(par.n),'_neqn_');
output_filename = strcat(output_filename,num2str(M),'.txt');
write_result(result,output_filename);

end

function f = ic(x,id)
if id == 1
    f = exp(-(x-0.5).*(x-0.5)*100);
    
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
            f = thetaIn * B(:,3); 
            
        case 2
            thetaIn = 0;
            f = thetaIn * B(:,3); 
    end

end


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

filename = strcat('result_Inflow2D/result_Reference/inflow_t_',num2str(t), ...
                '_points_',num2str(size(DX,2)),'_neqn_',num2str(length(U)),'.txt');

% norm_f means the norms of all the different moments
dlmwrite(filename,norm_f,'delimiter','\t','precision',10);

% norm_dx_f means the norm of the gradient of all the different moments
dlmwrite(filename,norm_dx_f,'delimiter','\t','-append','precision',10);
dlmwrite(filename,norm_dt_f,'delimiter','\t','-append','precision',10);
end
