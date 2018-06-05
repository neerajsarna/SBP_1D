function output = solver_DVM_3D(par)

if ~isfield(par,'save_during'), par.save_during = false; end % default value of save during the computation

if ~isfield(par,'ic'),       par.ic = @zero; end% Default: no init. cond.

if ~isfield(par,'source'),  par.source = @zero; end % Default: no source.

if par.num_bc ~=2
    assert(1 == 0, 'not valid num bc'); 
end   

% if we are not looking for the steady state then take the residual to be
% simply zero. 
if ~isfield(par,'steady_state')
    par.steady_state = false;
end

% find which of the given data is time dependent
time_dep = [nargin(par.source)>2 nargin(par.bc_inhomo)>7];
        
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

disp('delta t');
disp(par.dt);

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
bc_g = cell(par.num_bc,2);

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

% we have two systems now, one corresponding to g and the other to h
U = cell(2,par.n_eqn); dxU = U; UTemp = U; 
fM = cell(2);

fM{1} = zeros(par.n_eqn,length(X));
fM{2} = zeros(par.n_eqn,length(X));

% data structure for storing the values at the boundaries. The one coming
% from Sigma * B and the one coming from Sigma * g are both stored in
% bc_values. 
bc_values = cell(2,par.n_eqn);
k_RK = cell(4,2);

for i = 1 : 2
    for j = 1:par.n_eqn
        U{i,j} = X*0 + capargs(par.ic,X,j,par.Ax,par.Ay,i,par.mass_matrix);
        bc_values{i,j} = X*0;
    end
end

%% Time Loop
cputime = zeros(1,3);
t = 0; step_count = 0;

% compute the boundary inhomogeneity
for i = 1 : 2
    for j = 1:par.num_bc
        % need to convert to cell for the computations which follow
        bc_g{j,i} = num2cell(capargs(par.bc_inhomo,par.B{j},j,par.Ax, ...
                            par.Ay,i,U,par.all_w,t));
        
    end
end


residual = 0;

while t < par.t_end || residual > 10^(-5)
   
    residual = 0;
    
    if ~par.steady_state
        if t+par.dt > par.t_end
            par.dt = par.t_end-t;
        end
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
         
         
         if evaluate(2)
             for i = 1 : 2
                 for j = 1:par.num_bc
                     % need to convert to cell for the computations which follow
                     bc_g{j,i} = num2cell(capargs(par.bc_inhomo,par.B{j}, ...
                                            j,par.Ax,par.Ay,i,UTemp,par.all_w,t_temp(RK)));
                     
                 end
             end
         end
         
         rho = par.compute_density(UTemp,par.Ax,par.Ay,par.all_w);
         [ux,uy] = par.compute_velocity(UTemp,par.Ax,par.Ay,par.all_w);
         theta = par.compute_theta(UTemp,par.Ax,par.Ay,par.all_w);
        
         % we need to reconstruct the maxwellian
         for i = 1 : 2
             for k = 1 : length(X)
                 fM{i}(:,k) = minimize_entropy(par.Ax,par.Ay,...
                              par.mass_matrix,rho(k),ux(k),uy(k),theta(k),i);
             end
         end
             
         % compute the derivatives for g and h
         for i = 1 : 2
             for j = 1 : par.n_eqn
                 dxU{i,j} = DX * UTemp{i,j};
             end
             
          % extract all the value at x = x_start.
             bc_ID = 1;
             values = cellfun(@(a) a(1),UTemp(i,:),'Un',0);
             
             if par.pres_ID1
                 for j = 1 : par.n_eqn
                     % the term, values(bc_coupling{bc_ID}{j}), gives us the
                     % value of all the variables, at the boundary, which are
                     % coupled with the j-th variable.
                     % par.system_data.penalty_B{bc_ID}(j,bc_coupling{bc_ID}{j})
                     % gives us the j-th row of the penalty matrix and the
                     % entries in all those columns which have no zeros.
                     bc_values{i,j}(1) = bc_scaling * (sumcell(values(bc_coupling{bc_ID}{j}),...
                         par.penalty_B{bc_ID}(j,bc_coupling{bc_ID}{j})) - ...
                         sumcell(bc_g{bc_ID,i}(bc_coupling_g{bc_ID}{j}),...
                         par.penalty{bc_ID}(j,bc_coupling_g{bc_ID}{j})));
                 end
             end
             
             
             % extract all the value at x = x_end.
             bc_ID = 2;
             values = cellfun(@(a) a(end),UTemp(i,:),'Un',0);
             
             if par.pres_ID2
                 for j = 1 : par.n_eqn
                     bc_values{i,j}(end) = bc_scaling * ...
                         (sumcell(values(bc_coupling{bc_ID}{j}),...
                         par.penalty_B{bc_ID}(j,bc_coupling{bc_ID}{j})) - ...
                         sumcell(bc_g{bc_ID,i}(bc_coupling_g{bc_ID}{j}),...
                         par.penalty{bc_ID}(j,bc_coupling_g{bc_ID}{j})));
                 end
             end
             
            
             
             for j = 1 : par.n_eqn
                 
                 % multiplication of the derivatives and the system matrices
                 W = -sumcell(dxU(i,Ix{j}),par.Ax(j,Ix{j}));
                 
                 k_RK{RK,i}{j} = (W + bc_values{i,j});
             
                 k_RK{RK,i}{j} = k_RK{RK,i}{j} + ...
                                 (-UTemp{i,j}+fM{i}(j,:)')/par.Kn;
             end
         end
         
         % we should not update Utemp, we have not computed the two coupled
         % systems.
         for i = 1 : 2
             for j = 1 : par.n_eqn
                 if RK ~= 4
                     UTemp{i,j} = U{i,j} + k_RK{RK,i}{j} * dt_temp(RK + 1);
                     
                 end
             end
         end
     end
    
    for RK = 1 : 4
        for i = 1 : 2
            for j = 1 : par.n_eqn
                U{i,j} = U{i,j} + weight(RK) * k_RK{RK,i}{j} * par.dt;
                
                % if we are computing for the steady state
                if par.steady_state
                    residual = residual + norm(weight(RK) * k_RK{RK,i}{j})^2;
                end
            end
        end
    end
    
    % earlier we computed the square of the norm
    residual = sqrt(residual);
    
    
    
    step_count = step_count + 1;
    t = t + par.dt;
    cputime(1) = cputime(1) + toc;


    if mod(step_count,50) == 0
        disp('time: neqn: step_count: residual: ');
        disp(t);
        disp(par.n_eqn);
        disp(step_count);
        disp(residual);
    end
    
    tic
    
    if par.t_plot
        var_plot = par.compute_theta(U,par.Ax,par.Ay,par.all_w);
        %var_plot = cell2mat(cellfun(@(a) a(end),U,'Un',0));
        %plot(diag(par.Ax),var_plot,'-o');
        plot(X,var_plot,'-o');
        %xlim([min(par.x_m) max(par.x_p)]);
        xlim(par.ax);
        ylim([-1 1]);        
        drawnow;

%         v_id_plot = 1;
%         
%         % velocity corresponding to which we plot
%         v_plot = par.Ax(v_id_plot,v_id_plot);
%         
%         % the true solution shift
%         v_shift_plot = X - v_plot * t;
%         
%         exact_solution = zeros(length(v_shift_plot),1);
%         
%         for i = 1:length(v_shift_plot)
%             % take from the initial conditions
%             t_temp = (v_shift_plot(i)-1)/abs(v_plot);
%             if t_temp  <= 0
%                 exact_solution(i) = U{1,v_id_plot}(i);
%                 % take from the boundary conditions
%             else 
%                 % distance from the right boundary 
% %                 if t_temp <= 1
% %                     exact_solution(i) = exp(-1/(1-(t_temp-1)^2)) * exp(1);
% %                 else
% %                     exact_solution(i) = 1;
% %                 end
%                 
%                 exact_solution(i) = 0;
%             end
%         end
% 
%         
%         var_plot = U{1,v_id_plot};
%         plot(X,var_plot,'-o',X,exact_solution,'-*');
%         xlim(par.ax);
        
%        drawnow;
    end
    
    % option to save the data during time step, required for convergence
    % studies
    if par.save_during && mod(step_count,100) == 0
        % we compute the norms of the different features of the solution
        capargs(par.compute_during,U,weight,k_RK,PX,DX,t);
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

% we convert the data into contourf compatible format
function f = contourf_comp(Ax,X,data)
f = zeros(size(Ax,1),length(X));

for i = 1 : size(Ax,1)
    f(i,:) = data{i};
end

end

