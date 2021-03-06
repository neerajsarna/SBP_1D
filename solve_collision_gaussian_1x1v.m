function solve_collision_gaussian_1x1v(neqn)

par = struct(...
'name','1D Advection',... % name of example
'ic',@ic,... % initial conditions
'relax',@relax,... % production term defined below
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
'save_norms', @save_norms ...
);

par.t_plot = false;

par.n = 300;

par.n_eqn = neqn;

if par.n_eqn == 200
    par.save_during = true;
end
   
% we need the boundary matrix and the penalty matrix for both the
% boundaries
par.B = cell(par.num_bc,1);
par.penalty = cell(par.num_bc,1);

% first we develop the filenames
filenames_wall = dvlp_filenames_wall(par.n_eqn);
filenames_inflow = dvlp_filenames_inflow(par.n_eqn);

% we develop the system data
% par.Ax = dvlp_Ax1D(par.n_eqn);
% par.B{2} = dvlp_BInflow1D(par.n_eqn);
[par.Ax,Bwall] = get_system_data(filenames_wall);
[~,Binflow] = get_system_data(filenames_inflow);

% stabilise the boundary conditions, was slow in mathematica
Bwall = stabilize_boundary(par.Ax,Bwall);
Binflow = stabilize_boundary(par.Ax,Binflow);

par.B{2} = Bwall;

% develop the boundary matrix at ID2 or x=0
par.B{1} = dvlp_B_ID1(Bwall);

% develop the characteristic penalty matrix
[par.penalty{1}] = dvlp_penalty(-par.Ax,par.B{1});
[par.penalty{2}] = dvlp_penalty(par.Ax,par.B{2});

% develop the penalty*B matrix
for i = 1 : par.num_bc
    par.penalty_B{i} = par.penalty{i} * par.B{i};
end

result = solver(par);

output_filename = strcat('result_collision_gaussian_1x1v/wall_tend_',...
                        num2str(par.t_end),'_points_',num2str(par.n),'_neqn_');
disp('written to:');
disp(output_filename);
output_filename = strcat(output_filename,num2str(par.n_eqn),'.txt');
write_result(result,output_filename);

end


% working on inflow boundaries, we consider vacum boundary conditions
function f = bc_inhomo(B,bc_id,t)

    switch bc_id
        % boundary at x = x_start
        case 1
            thetaIn = 0;
           
            f = thetaIn * B(:,3)/sqrt(2); 
            
         % x = 1
        case 2
            thetaIn = 0;
            f = thetaIn * B(:,3)/sqrt(2); 
    end

end
        
function f = ic(x,id)

density = exp(-(x-0.5).*(x-0.5)*100);
velocity = exp(-(x-0.5).*(x-0.5)*100);

if id == 1
    f = density;
    % f = density * 0;
end
if id == 2
    f = velocity;
end

if id > 2
    f = density * 0;
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

function[filenames] = dvlp_filenames_inflow(nEqn)

filenames = struct;

filenames.B = strcat("system_matrices1D/Binflow_1D_",num2str(nEqn));
filenames.B = strcat(filenames.B,".txt");

filenames.Ax = strcat("system_matrices1D/A1_1D_",num2str(nEqn));
filenames.Ax = strcat(filenames.Ax,".txt");
end

function[filenames] = dvlp_filenames_wall(nEqn)

filenames = struct;

filenames.B = strcat("system_matrices1D/Bwall_1D_",num2str(nEqn));
filenames.B = strcat(filenames.B,".txt");

filenames.Ax = strcat("system_matrices1D/A1_1D_",num2str(nEqn));
filenames.Ax = strcat(filenames.Ax,".txt");
end


function [int_f,int_dx_f,int_dt_f] = compute_during(U,weight,k_RK,PX,DX,t,t_Old)
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

% we perform the integral here itself 
% using the midpoint rule
int_f = norm_f * (t-t_Old);
int_dx_f = norm_dx_f * (t-t_Old);
int_dt_f = norm_dt_f * (t-t_Old);
end

function save_norms(int_f,int_dx_f,int_dt_f,n)
filename = strcat('result_collision_gaussian_1x1v/result_Reference/wall_norms', ...
                '_points_',num2str(n),'_neqn_',num2str(length(int_f)),'.txt');

dlmwrite(filename,int_f,'delimiter','\t','precision',10);
dlmwrite(filename,int_dx_f,'delimiter','\t','-append','precision',10);
dlmwrite(filename,int_dt_f,'delimiter','\t','-append','precision',10);
end
