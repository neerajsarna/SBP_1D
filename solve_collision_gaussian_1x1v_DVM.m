% we solve the inflow problem using discrete velocity method 
function solve_collision_gaussian_1x1v_DVM(nc)

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

par.n = 300;
[par.x_m,par.w_m] = gauss_quadrature(nc,-5,0);
[par.x_p,par.w_p] = gauss_quadrature(nc,0,5);

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

result= solver_DVM_1x1v(par);

temp = cell(length(result),1);

for i = 1 : length(temp)
    temp{i} = result(i).sol;
end

density = compute_density(temp,par.Ax,par.all_w);
u = compute_velocity(temp,par.Ax,par.all_w);
alpha2 = compute_alpha2(temp,par.Ax,par.all_w);

filename = strcat('result_collision_gaussian_1x1v/DVM_wall_tend_',num2str(par.t_end),...
                        '_points_',num2str(par.n),'_neqn_');
filename = strcat(filename,num2str(nc),'.txt');

dlmwrite(filename,result(1).X','delimiter','\t','precision',10);
dlmwrite(filename,density','delimiter','\t','-append','precision',10);
dlmwrite(filename,u','delimiter','\t','-append','precision',10);
dlmwrite(filename,alpha2','delimiter','\t','-append','precision',10);

% a plot of the distribution function along the x and the velocity space
temp = zeros(length(result),length(result(1).X));

[v_grid,sort_grid] = sort([par.x_m',par.x_p']);

for i = 1:length(result)
    temp(i,:) = result(sort_grid(i)).sol;
end

[x_mesh,v_mesh] = meshgrid(result(1).X,v_grid);

%surf(x_mesh,v_mesh,temp);

output_filename = strcat('result_collision_gaussian_1x1v/DVM_f_wall_tend_', ...
                        num2str(par.t_end),'_points_',num2str(par.n),'_neqn_');
                    
output_filename = strcat(output_filename,num2str(nc),'.txt');

dlmwrite(output_filename,result(1).X','delimiter','\t','precision',10);
dlmwrite(output_filename,v_grid,'-append','delimiter','\t','precision',10);
dlmwrite(output_filename,par.all_w(sort_grid),'-append','delimiter','\t','precision',10);
dlmwrite(output_filename,temp,'-append','delimiter','\t','precision',10);
% 
% output_filename = strcat('result_Inflow_fluctuateT_1x1v/PX.txt');
% dlmwrite(output_filename,full(result(1).P),'delimiter','\t','precision',10);
% 
% output_filename = strcat('result_Inflow_fluctuateT_1x1v/DX.txt');
% dlmwrite(output_filename,full(result(1).D),'delimiter','\t','precision',10);

end


% working on inflow boundaries, we consider vacum boundary conditions
% the inhomogeneity is time dependent because it depends upon U
function f = bc_inhomo(B,bc_id,Ax,U,all_weights,t)

    uW = 0;
    thetaW = 0;
    alpha2W = thetaW/sqrt(2);
    id = find(diag(B) == 1);
    f = B(:,1) * 0;
    
    switch bc_id
        % boundary at x = x_start
        case 1
            temp = cell2mat(cellfun(@(a) a(1),U,'Un',0));
            rhoW = compute_rhoW(Ax,temp,thetaW,all_weights,bc_id);        
        case 2
            temp = cell2mat(cellfun(@(a) a(end),U,'Un',0));
            rhoW = compute_rhoW(Ax,temp,thetaW,all_weights,bc_id);
    end
    
    for i = 1 : length(id)
        index = id(i);
        f(index) = compute_fM(Ax,rhoW,uW,alpha2W,index);
    end

end
        
function f = ic(x,id,Ax)
density = exp(-(x-0.5).*(x-0.5)*100);
velocity = exp(-(x-0.5).*(x-0.5)*100);

% the velocity corresponding to the id
v = Ax(id,id);

%f = density * f0(v);
f = f0(v) * (density + v * velocity);
end

function f = f0(v)
f = exp(-v.^2/2)/sqrt(2 * pi);
end

function f = relax(x,id)
Kn = 0.1;
f = zeros(length(x),1);

end

function f = compute_density(U,Ax,all_weights)
grid_points = length(U{1});

f = zeros(grid_points,1);

% loop over the grid
for j = 1 : length(all_weights)
    for i = 1 : grid_points
        f(i) = all_weights(j) * U{j}(i) + f(i);
    end
end

end

function f = compute_velocity(U,Ax,all_weights)
grid_points = length(U{1});

f = zeros(grid_points,1);

% scale by velocity
all_weights = diag(Ax)'.*all_weights;

% loop over the grid
for j = 1 : length(all_weights)
    for i = 1 : grid_points
        f(i) = all_weights(j) * U{j}(i) + f(i);
    end
end

end

% second coefficient in the expansion
function f = compute_alpha2(U,Ax,all_weights)
grid_points = length(U{1});

f = zeros(grid_points,1);

% scale by He_2
all_weights = (diag(Ax)'.*diag(Ax)'-1).*all_weights/sqrt(2);

% loop over the grid
for j = 1 : length(all_weights)
    for i = 1 : grid_points
        f(i) = all_weights(j) * U{j}(i) + f(i);
    end
end

end

% compute the Maxwellian 
function f = compute_fM(Ax,rho,u,alpha2,id)
v = Ax(id,id);

f = f0(v) * (rho + u * v + alpha2 * (v^2-1)/sqrt(2));
end

function rhoW = compute_rhoW(Ax,U,thetaW,all_weights,bc_id)
U = U';

all_vx = diag(Ax);
pos_vx_p = find(all_vx > 0);
pos_vx_m = find(all_vx < 0);

vx_p = all_vx(pos_vx_p);
vx_m = all_vx(pos_vx_m);

num_pos = length(vx_p);
num_neg = length(vx_m);

weight_vx_m = all_weights(pos_vx_m);
weight_vx_p = all_weights(pos_vx_p);

int_f0 = 0;
int_f0_sq_vx = 0;
int_f_vx = 0;

switch bc_id
    
    % x = 1
    case 2
        % integrate the maxwellian
        for i = 1 : num_neg
            % half integral of f0 * vx
            int_f0 = int_f0 + vx_m(i) * f0(vx_m(i)) * weight_vx_m(i);
            % half integral of f0 * vx^2 * vx
            int_f0_sq_vx = int_f0_sq_vx...
                          + vx_m(i) * vx_m(i)^2 * f0(vx_m(i)) * weight_vx_m(i);
        end
        
        % integrate the kinetic solution
        for i = 1 : num_pos
            int_f_vx = vx_p(i)*U(pos_vx_p(i))*weight_vx_p(i) + int_f_vx;
        end
       
    case 1
        % integrate the maxwellian
        for i = 1 : num_pos
            int_f0 = int_f0 + vx_p(i) * f0(vx_p(i)) * weight_vx_p(i);
            int_f0_sq_vx = int_f0_sq_vx ...
                        + vx_p(i) * vx_p(i)^2 * f0(vx_p(i)) * weight_vx_p(i);
        end
        
        for i = 1 : num_neg
            int_f_vx = vx_m(i)*U(pos_vx_m(i))*weight_vx_m(i) + int_f_vx;
        end
        
end
        rhoW = -int_f_vx-thetaW * (int_f0_sq_vx-int_f0)/2;
        rhoW = rhoW/int_f0;
end
