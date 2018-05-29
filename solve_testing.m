function solve_testing(neqn)

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
'exact_solution',@exact_solution...
);

par.t_plot = true;

par.n = 300;

par.n_eqn = neqn;

   
% toy advection equation
par.Ax = -2 * eye(neqn);

% at x = 1 we prescribe boundary to negative velocities
par.B{2} = diag(double(diag(par.Ax)<0));
par.B{1} = diag(double(diag(-par.Ax)<0));

% prescribe a value to the negative velocities
par.penalty{2} = (par.Ax-abs(par.Ax))/2;

% characteristic penalty matrix, the matrix X_n^{-} will just be identity
par.penalty{1} = (-par.Ax-abs(par.Ax))/2;

% develop the penalty*B matrix
for i = 1 : par.num_bc
    par.penalty_B{i} = par.penalty{i} * par.B{i};
end

result = solver(par);

end


% working on inflow boundaries, we consider vacum boundary conditions
function f = bc_inhomo(B,bc_id,t)
    if t <= 1
        thetaIn = exp(-1/(1-(t-1)^2)) * exp(1);
    else
        thetaIn = 1;
    end
  
    
    f = diag(B) * thetaIn;

end
        
function f = ic(x,id)
if id == 1
    f = exp(-(x-0.5).*(x-0.5)*100);
    
else
    f = zeros(length(x),1);
end
end

function f = relax(x,id)
Kn = inf;
f = zeros(length(x),1);

% anything above temperature has to be relaxed
if id > 3
    f = -1/Kn * ones(length(x),1);
end

end

function f = exact_solution(X,t,Ax,U)
        
        v_plot = Ax(1,1);
        
        % the true solution shift
        v_shift_plot = X - v_plot * t;
        
        f = zeros(length(v_shift_plot),1);
        
        for i = 1:length(v_shift_plot)
            % take from the initial conditions
            t_temp = (v_shift_plot(i)-1)/2;
            if t_temp  <= 0
                f(i) = U{1}(i);
                % take from the boundary conditions
            else 
                % distance from the right boundary 
                if t_temp <= 1
                    f(i) = exp(-1/(1-(t_temp-1)^2)) * exp(1);
                else
                    f(i) = 1;
                end
            end
        end
end




