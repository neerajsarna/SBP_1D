clear all;

par = struct(...
'name','1D Advection',... % name of example
'n_eqn',1,... % number of variables in the system
'source',@source,... % source term (defined below)
'ic',@ic,... % initial conditions
'bc_inhomo',@bc_inhomo,... % source term (defined below)
'ax',[0 1],... % coordinates of computational domain
 't_end',1.0,... % the end time of the computatio
 'CFL',2.0,...      % the crude cfl number
 'num_bc',2,... % number of boundaries in the domain
 'num_bc_pres',1,... % number of boundaries conditions to be prescribed
 'pres_ID1',true,... % whether we need to prescribe something at x = x_start
 'pres_ID2',false,... % whether we need to prescribe something at x = x_start
 'var_output',1,... % the variable which should be output
'output',@output... % problem-specific output routine (defined below
);

par.t_plot = false;

par.Ax = 1;

% the boundary matrix
par.B{1} = 1;
par.B{2} = 1;

% at ID = 1, we prescribe the boundary conditions
par.penalty_B{1} = -1;
par.penalty{1} = -1;

% at ID = 2, nothing happens
par.penalty_B{2} = 0;
par.penalty{2} = 0;


error = [];
grid_spacing = [];
for i = 2:8
    
    par.n = 2^i;
    
    disp('solving for :');
    disp(par.n);
    
    result = solver(par);
    
    error = [error sqrt((result.U-exact_solution(result.X,par.t_end))' * result.P * ...
                (result.U-exact_solution(result.X,par.t_end)))];
            
    grid_spacing = [grid_spacing result.h];    
        
end

convg_order = log(error(end)/error(1))/log(grid_spacing(end)/grid_spacing(1));
disp('convergence order');
disp(convg_order);

function f = source(x,index)
    f = 0;
end

function f = ic(x,id)
    switch id
        case 1
            f = sin(pi * x);
    end
end

function f = bc_inhomo(B,bc_id,t)

switch bc_id
    case 1
    f = -sin(pi*t);
    case 2
    f = 0;
end
end

function f = exact_solution(x,t)
f = sin(pi * (x-t));
end
