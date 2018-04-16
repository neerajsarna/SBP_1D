function output = solver(par)

if ~isfield(par,'save_during'), par.save_during = false; end % default value of save during the computation

if ~isfield(par,'ic'),       par.ic = @zero; end% Default: no init. cond.


if ~isfield(par,'relax'),       par.relax = @zero; end% Default: no init. cond.


if ~isfield(par,'source'),  par.source = @zero; end % Default: no source.


if par.num_bc ~=2
    assert(1 == 0, 'not valid num bc'); 
end   

% find which of the given data is time dependent
time_dep = [nargin(par.source)>2 nargin(par.bc_inhomo)>2];
        
% if the number of arguments are greater than 3 then definitely we have 
% an anisotropic source term.
if ~isfield(par,'source_ind') 
    par.source_ind = 1:par.n_eqn;  
end

% corresponding to every row in Ax, stores the non-zero indices
Ix = cellfun(@find,num2cell(par.Ax',1),'Un',0);         

% size of the grid
h = (par.ax(2)-par.ax(1))/par.n;

% a crude approximation for delta_t
par.dt = min(h)/abs(eigs(par.Ax,1,'lm'))/par.CFL;

% the grid in the x-direction 
X = par.ax(1):h(1):par.ax(2);
X = X';

% we need to know which elemtns are coupled with which one at the
% boundaries. The 
% ID = 1
% a loop over all the boundaries
% consider the coupling for the Sigma * B term
bc_coupling = cell(par.num_bc,1);

% consider the coupling for the Sigma  * g term
bc_coupling_g = cell(par.num_bc,1);

%the boundary inhomogeneity. Every component of the cell will be equal to
%the number of boundar conditions which we need to prescribe. 
bc_g = cell(par.num_bc,1);

% we need to do this fixing for the machine error
for i = 1 : par.num_bc
    % for the boundary id=i, we find the variables to which every variable
    % in the moment system
    % is coupled 
    % the matrix penalty_B{i} contains Sigma * B at the boundary with ID =
    % i. the matrix penalty{i} contains the penalty matrix at the boundary
    % with ID = i.
    bc_coupling{i} = cellfun(@(a) find(abs(a) > 1e-14) ,num2cell(par.penalty_B{i},2),'Un',0);
    bc_coupling_g{i} = cellfun(@(a) find(abs(a)> 1e-14) ,num2cell(par.penalty{i},2),'Un',0);
    

end

[DX,PX,InvPX] = sbp_traditional_2(h,par.n);

% scaling for the boundary conditions. See the SBP paper. 
bc_scaling = InvPX(1,1);

% initialize the solution variables
% n_sys is the total number of equations in the system
% dxU and dyU are the derivatives in the x and y direction.
% force is the forcing term 
% set the initial conditions in case they are needed else initialize with
% zero
U = cell(1,par.n_eqn); dxU = U; force = U; UTemp = U; Prod = U;

% data structure for storing the values at the boundaries. The one coming
% from Sigma * B and the one coming from Sigma * g are both stored in
% bc_values. 
bc_values = cell(1,par.n_eqn);
k_RK = cell(4,1);
for j = 1:par.n_eqn                          
    U{j} = X*0 + capargs(par.ic,X,j);     
    force{j} = X*0;
    bc_values{j} = X*0;
    Prod{j} = X*0;
end

%% Time Loop
cputime = zeros(1,3);
t = 0; step_count = 0;


for j = par.source_ind     
        % apply source to the active components. We store the value on the
        % whole grid at once. 
        force{j} = capargs(par.source,X,j,t); 
end

% compute the boundary inhomogeneity
for j = 1:par.num_bc
    % need to convert to cell for the computations which follow
    bc_g{j} = num2cell(capargs(par.bc_inhomo,par.B{j},j,t));
   
end

for j = 1 : par.n_eqn
    Prod{j} = capargs(par.relax,X,j);
end


while t < par.t_end
    
    if t+par.dt > par.t_end
        par.dt = par.t_end-t;
    end
    
    % the ode sytem can be written as U_t = Op.
    % RK = 4 implementation
     tic
            t_temp = [t t + par.dt/2 t + par.dt/2 t + par.dt];
            dt_temp = [0 par.dt/2 par.dt/2 par.dt];
            weight = [1/6 2/6 2/6 1/6];
     
     UTemp = U; 
     for RK = 1 : 4
            evaluate = time_dep & (t_temp(RK) > 0);
            
            if evaluate(1)
                for j = par.source_ind    
                         force{j} = capargs(par.source,X,j,t_temp(RK)); 
                end
            end
            
            if evaluate(2)
                for j = 1:par.num_bc
                 % need to convert to cell for the computations which follow
                    bc_g{j} = num2cell(capargs(par.bc_inhomo,par.B{j},j,t_temp(RK)));
   
                end
            end
            
            
            for i = 1:par.n_eqn
                dxU{i} = DX * UTemp{i};
            end

            
            % extract all the value at x = x_start. 
            bc_ID = 1;
            values = cellfun(@(a) a(1),UTemp,'Un',0);
            
            if par.pres_ID1
            for j = 1 : par.n_eqn
%                 the term, values(bc_coupling{bc_ID}{j}), gives us the
%                 value of all the variables, at the boundary, which are
%                 coupled with the j-th variable. 
%                 par.system_data.penalty_B{bc_ID}(j,bc_coupling{bc_ID}{j})
%                 gives us the j-th row of the penalty matrix and the
%                 entries in all those columns which have no zeros. 
                bc_values{j}(1) = bc_scaling * (sumcell(values(bc_coupling{bc_ID}{j}),...
                                      par.penalty_B{bc_ID}(j,bc_coupling{bc_ID}{j})) - ...
                                      sumcell(bc_g{bc_ID}(bc_coupling_g{bc_ID}{j}),...
                                      par.penalty{bc_ID}(j,bc_coupling_g{bc_ID}{j})));
            end
            end
            

           % extract all the value at x = x_end. 
            bc_ID = 2;
            values = cellfun(@(a) a(end),UTemp,'Un',0);
            
            if par.pres_ID2
            for j = 1 : par.n_eqn    
                bc_values{j}(end) = bc_scaling * (sumcell(values(bc_coupling{bc_ID}{j}),...
                                                 par.penalty_B{bc_ID}(j,bc_coupling{bc_ID}{j})) - ...
                                                 sumcell(bc_g{bc_ID}(bc_coupling_g{bc_ID}{j}),...
                                                 par.penalty{bc_ID}(j,bc_coupling_g{bc_ID}{j})));
            end
            end
            
            % checking conservation 
%             flux_wall = ones(1,length(X)) * PX * (bc_values{2});
%             flux_int = -ones(1,length(X)) * PX * (dxU{1} + dxU{3} * sqrt(2));
%             
%             disp('velocity norm');
%             disp(norm(UTemp{2},2));
            
            for i = 1 : par.n_eqn
                
                % multiplication of the derivatives and the system matrices
                W = -sumcell(dxU(Ix{i}),par.Ax(i,Ix{i}));
                                    
                k_RK{RK}{i} = (W + force{i} + bc_values{i} + Prod{i}.*UTemp{i});
                
                if RK ~= 4
                    UTemp{i} = U{i} + k_RK{RK}{i} * dt_temp(RK + 1);
                end
            end
                
     end
    
    for RK = 1 : 4
        for i = 1 : par.n_eqn
            U{i} = U{i} + weight(RK) * k_RK{RK}{i} * par.dt;
        end
    end
    
    step_count = step_count + 1;
    t = t + par.dt;
    cputime(1) = cputime(1) + toc;
    
    if mod(step_count,50) == 0
        disp('time: neqn: ');
        disp(t);
        disp(par.n_eqn);
        disp(step_count);
    end
    
    tic
    
    if par.t_plot
        plot(X,U{par.var_output},'-o');
        xlim(par.ax);
        ylim([-1 1]);        
        drawnow;
    end
    
    % option to save the data during time step, required for convergence
    % studies
    if par.save_during && mod(step_count,500) == 0
        output = struct('X',X, ...
                        'sol',U, ...
                        'P',PX, ...
                        'h',h);

        output_filename = strcat('result_wall_specular/result_Reference/inflow_t_',num2str(t),'_points_',num2str(par.n),'_neqn_');
        output_filename = strcat(output_filename,num2str(par.n_eqn),'.txt');
        write_result(output,output_filename);
    end
    
    cputime(2) = cputime(2) + toc;
end

fprintf('%0.0f time steps\n',step_count)           % Display test
cputime = reshape([cputime;cputime/sum(cputime)*1e2],1,[]);   % case info
fprintf(['CPU-times\n advection:%15.2fs%5.0f%%\n',... % and CPU times.
    'plotting:%16.2fs%5.0f%%\n'],cputime)

output = struct('X',X, ...
                'sol',U, ...
                 'P',PX, ...
                 'h',h);
end


function z = capargs(fct,varargin)
% Call function fct with as many arguments as it requires (at least 1),
% and ignore further arguments.
narg = max(nargin(fct),1);
z = fct(varargin{1:narg});

end

function f = zero(varargin)
% Zero function.
f = zeros(size(varargin{1}));
end

function S = sumcell(A,w)
% Add vector of cells A, weighted with vector w.
S = 0;

for j = 1:length(w)
    S = S+A{j}*w(j);
end

end

