% manufactured solution for the 1D moment system
clear all;

par = struct(...
'name','moments1D Manufactured Solutions',... % name of example
'ic',@ic,... % initial conditions
'relax',@relax,... % production term defined below
'source',@source,... % defined below
'bc_inhomo',@bc_inhomo,... % source term (defined below)
'ax',[0 1],... % coordinates of computational domain
 't_end',0.5,... % the end time of the computatio
 'CFL',4.0,...      % the crude cfl number
 'num_bc',2,... % number of boundaries in the domain
 'pres_ID1',true,... % whether we need to prescribe something at x = x_start
 'pres_ID2',true,... % whether we need to prescribe something at x = x_end
 'var_output',1,... % the variable which should be plotted
'output',@output... % problem-specific output routine (defined below)
);

par.t_plot = false;

par.n_eqn = 4;

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

errorTot = [];
gridSpace = [];

for i = 3:8
    par.n = 2^i;
    disp('solving for :');
    disp(2^i);
    result = solver(par);
    
    error = 0;
    
    for j = 1 : length(result)
        error = error + (result(j).sol-exact_solution(result(j).X,j,par.t_end))'*result(j).P * ...
                (result(j).sol-exact_solution(result(j).X,j,par.t_end));   
    end
    
    errorTot = [errorTot sqrt(error)];
    gridSpace = [gridSpace result(1).h];
end

expect_order = 2.0;
ref_Order = exact_order(gridSpace(1),errorTot(1),gridSpace(end),expect_order);
loglog(gridSpace,errorTot,'-o',ref_Order([1 2]),ref_Order([3 4]),'->');
grid on;
disp('convg order');
disp(log(errorTot(end)/errorTot(1))/log(gridSpace(end)/gridSpace(1)));
    
function f = bc_inhomo(B,bc_id,t)
    
neqn = size(B,2);

    switch bc_id

        % compute the exact solution at x = 1
        case 1            
            exact_sol = arrayfun(@exact_solution,0 * ones(1,neqn),1:neqn,t*ones(1,neqn));
            
            
        case 2
            exact_sol = arrayfun(@exact_solution,1*ones(1,neqn),1:neqn,t*ones(1,neqn));
    end
    
    f = B*exact_sol'; 

end
        
function f = ic(x,id)

f = sin(pi * x/id);

end

function f = relax(x,id)
Kn = 1.0;
f = zeros(length(x),1);

% anything above temperature has to be relaxed
if id > 3
    f = -1/Kn * ones(length(x),1);
end

end

% the manufactured solution
function f = exact_solution(x,id,t)
    f = sin(pi * x/id) * cos(pi * t/id);
end

% the source term
function f = source(x,id,t)
    switch id
        case 1
            f = (pi*cos((pi*t)/2.)*cos((pi*x)/2.))/2. - pi*sin(pi*t)*sin(pi*x);
        case 2
            f = (sqrt(2)*pi*cos((pi*t)/3.)*cos((pi*x)/3.))/3. + pi*cos(pi*t)*cos(pi*x) - ...
                (pi*sin((pi*t)/2.)*sin((pi*x)/2.))/2.;
        case 3
            f = (sqrt(3)*pi*cos((pi*t)/4.)*cos((pi*x)/4.))/4. + (pi*cos((pi*t)/2.)*cos((pi*x)/2.))/sqrt(2) - ...
                (pi*sin((pi*t)/3.)*sin((pi*x)/3.))/3.;
        case 4
            f = (pi*cos((pi*t)/3.)*cos((pi*x)/3.))/sqrt(3) + cos((pi*t)/4.)*sin((pi*x)/4.) - ...
                (pi*sin((pi*t)/4.)*sin((pi*x)/4.))/4.;
    end
end

function[filenames] = dvlp_filenames(nEqn)

filenames = struct;

filenames.B = strcat("generic_1D/system_matrices1D/Bwall_1D_",num2str(nEqn));
filenames.B = strcat(filenames.B,".txt");

filenames.Ax = strcat("generic_1D/system_matrices1D/A1_1D_",num2str(nEqn));
filenames.Ax = strcat(filenames.Ax,".txt");
end


