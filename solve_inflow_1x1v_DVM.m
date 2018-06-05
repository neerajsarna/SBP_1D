% we solve the inflow problem using discrete velocity method 
function solve_inflow_DVM(nc)

par = struct(...
'name','Inflow DVM',... % name of example
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
'compute_density',@compute_density, ...
'compute_velocity',@compute_velocity, ...
'compute_alpha2',@compute_alpha2, ...
'compute_fM',@compute_fM ...
);

par.Kn = 0.1;
par.t_plot = false;
par.n_eqn = 2 * nc;
%par.n_eqn =nc;

par.n = 100;
%[temp_x,temp_w] = gauss_quadrature(nc,-3,3);
[par.x_m,par.w_m] = gauss_quadrature(nc,-3,0);
[par.x_p,par.w_p] = gauss_quadrature(nc,0,3);

% par.x_m = temp_x(temp_x < 0);
% par.x_p = temp_x(temp_x >= 0);
% par.w_m = temp_w(temp_x < 0);
% par.w_p = temp_w(temp_x >= 0);

par.Ax = diag([par.x_m',par.x_p']);
par.all_w = [par.w_m',par.w_p'];

par.B{2} = diag([ones(1,length(par.x_m)),zeros(1,length(par.x_p))]);
par.B{1} = diag([zeros(1,length(par.x_m)),ones(1,length(par.x_p))]);

% prescribe a value to the negative velocities
par.penalty{2} = diag([par.x_m',zeros(1,length(par.x_p))]);

% prescribe a value to the positive velocities
par.penalty{1} = diag([zeros(1,length(par.x_m)),-par.x_p']);

% develop the penalty*B matrix
for i = 1 : par.num_bc
    par.penalty_B{i} = par.penalty{i} * par.B{i};
end

result= solver_DVM(par);

% temp = cell(length(result),1);
% 
% for i = 1 : length(temp)
%     temp{i} = result(i).sol;
% end
% 
% density = compute_density(temp,par.Ax,par.all_w);
% u = compute_velocity(temp,par.Ax,par.all_w);
% alpha2 = compute_alpha2(temp,par.Ax,par.all_w);

% filename = strcat('result_Inflow/DVM_inflow_tend_',num2str(par.t_end),...
%                         '_points_',num2str(par.n),'_neqn_');
% filename = strcat(filename,num2str(nc),'.txt');
% 
% dlmwrite(filename,result(1).X','delimiter','\t','precision',10);
% dlmwrite(filename,density','delimiter','\t','-append','precision',10);
% dlmwrite(filename,u','delimiter','\t','-append','precision',10);
% dlmwrite(filename,alpha2','delimiter','\t','-append','precision',10);

% a plot of the distribution function along the x and the velocity space
temp = zeros(length(result),length(result(1).X));

[v_grid,sort_grid] = sort([par.x_m',par.x_p']);

for i = 1:length(result)
    temp(i,:) = result(sort_grid(i)).sol;
end

[x_mesh,v_mesh] = meshgrid(result(1).X,v_grid);

surf(x_mesh,v_mesh,temp);

output_filename = strcat('result_Inflow/DVM_inflow_', ...
                        num2str(par.t_end),'_points_',num2str(par.n),'_neqn_');
output_filename = strcat(output_filename,num2str(nc),'.txt');

dlmwrite(output_filename,result(1).X','delimiter','\t','precision',10);
dlmwrite(output_filename,v_grid,'-append','delimiter','\t','precision',10);
dlmwrite(output_filename,temp,'-append','delimiter','\t','precision',10);

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
        
function f = ic(x,id,Ax)
density = exp(-(x-0.5).*(x-0.5)*100);

% the velocity corresponding to the id
v = Ax(id,id);

f = density * f0(v);
end

function f = f0(v)
f = exp(-v^2/2)/sqrt(2 * pi);
end

function f = relax(x,id)
Kn = 0.1;
f = zeros(length(x),1);

end

function f = compute_density(U,Ax,all_weights)
grid_points = length(U{1});

f = zeros(grid_points,1);

% loop over the grid
for i = 1 : grid_points
    temp = cell2mat(cellfun(@(a) a(i),U,'Un',0));
    % might need to reshape
    if size(temp,1) ~= size(all_weights,1)
        temp = reshape(temp,size(all_weights));
    end
    f(i) = sum(all_weights.*temp);
end

end

function f = compute_velocity(U,Ax,all_weights)
grid_points = length(U{1});

f = zeros(grid_points,1);

% scale by velocity
all_weights = diag(Ax)'.*all_weights;

% loop over the grid
for i = 1 : grid_points
    temp = cell2mat(cellfun(@(a) a(i),U,'Un',0));
    if size(temp,1) ~= size(all_weights,1)
        temp = reshape(temp,size(all_weights));
    end
    f(i) = sum(all_weights.*temp);
end

end

% second coefficient in the expansion
function f = compute_alpha2(U,Ax,all_weights)
grid_points = length(U{1});

f = zeros(grid_points,1);

% scale by velocity
all_weights = (diag(Ax)'.*diag(Ax)'-1).*all_weights/sqrt(2);

% loop over the grid
for i = 1 : grid_points
    temp = cell2mat(cellfun(@(a) a(i),U,'Un',0));
    if size(temp,1) ~= size(all_weights,1)
        temp = reshape(temp,size(all_weights));
    end
    f(i) = sum(all_weights.*temp);
end

end

% compute the Maxwellian 
function f = compute_fM(Ax,rho,u,alpha2,id)
v = Ax(id,id);

f = f0(v) * (rho + u * v + alpha2 * (v^2-1)/sqrt(2));
end
