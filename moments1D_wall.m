clear all;

par = struct(...
'name','1D Advection',... % name of example
'ic',@ic,... % initial conditions
'relax',@relax,... % production term defined below
'bc_inhomo',@bc_inhomo,... % source term (defined below)
'ax',[0 1],... % coordinates of computational domain
 't_end',10,... % the end time of the computatio
 'CFL',4.0,...      % the crude cfl number
 'num_bc',2,... % number of boundaries in the domain
 'pres_ID1',true,... % whether we need to prescribe something at x = x_start
 'pres_ID2',true,... % whether we need to prescribe something at x = x_end
 'var_output',3,... % the variable which should be plotted
'output',@output... % problem-specific output routine (defined below)
);

par.t_plot = false;

par.n_eqn = 4;

par.n = 50;

% we need the boundary matrix and the penalty matrix for both the
% boundaries
par.B = cell(par.num_bc,1);
par.penalty = cell(par.num_bc,1);

% first we develop the filenames
filenames = dvlp_filenames(par.n_eqn);

% we develop the system data
[par.Ax,par.B{2}] = get_system_data(filenames);

% develop the boundary matrix at ID2 or x=0
par.B{1} = dvlp_B_ID1(par.B{2});

% develop the characteristic penalty matrix
[par.penalty{1}] = dvlp_penalty(-par.Ax,par.B{1});
[par.penalty{2}] = dvlp_penalty(par.Ax,par.B{2});

% develop the penalty*B matrix
for i = 1 : par.num_bc
    par.penalty_B{i} = par.penalty{i} * par.B{i};
end

result = solver(par);

function f = bc_inhomo(B,bc_id)
    switch bc_id
        % boundary at x = x_start
        case 1
            thetaIn = -1;
            f = thetaIn * B(:,3)/sqrt(2); 
            
        case 2
            thetaIn = 1;
            f = thetaIn * B(:,3)/sqrt(2); 
    end

end
        
function f = ic(x,id)
% if id == 1
%     f = exp(-(x-0.5).^2 * 100);
% else
%     f = zeros(1,length(x));
% end

f = zeros(length(x),1);

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

filenames.B = strcat("generic_1D/system_matrices1D/Bwall_1D_",num2str(nEqn));
filenames.B = strcat(filenames.B,".txt");

filenames.Ax = strcat("generic_1D/system_matrices1D/A1_1D_",num2str(nEqn));
filenames.Ax = strcat(filenames.Ax,".txt");
end


