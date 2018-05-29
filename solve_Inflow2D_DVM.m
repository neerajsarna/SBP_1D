% we solve the Heat Conduction problem using discrete velocity method 
function solve_Inflow2D_DVM(nc)

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
'compute_theta',@compute_theta, ...
'compute_fM',@compute_fM ...
);

par.Kn = 0.1;
par.t_plot = false;
par.n_eqn = (2 * nc) * (2 * nc);
%par.n_eqn =nc;

% number of points in the spatial discretization
par.n = 50;

%[temp_x,temp_w] = gauss_quadrature(nc,-3,3);
[par.x_m,par.w_m] = gauss_quadrature(nc,-5,0);
[par.x_p,par.w_p] = gauss_quadrature(nc,0,5);

% velocity grid in one dimension
[v_grid1D,perm] = sort([par.x_m' par.x_p']);
w_grid1D = [par.w_m',par.w_p'];
w_grid1D = w_grid1D(perm);

% velocity grid in two dimensions
[vx_grid2D,vy_grid2D] = meshgrid(v_grid1D,v_grid1D);

% the quadrature weights grid in two dimensions
[wx_grid2D,wy_grid2D] = meshgrid(w_grid1D,w_grid1D);

w = wx_grid2D.*wy_grid2D;

% par.x_m = temp_x(temp_x < 0);
% par.x_p = temp_x(temp_x >= 0);
% par.w_m = temp_w(temp_x < 0);
% par.w_p = temp_w(temp_x >= 0);

par.Ax = diag(vx_grid2D(:));
par.Ay = diag(vy_grid2D(:));
par.all_w = w(:);


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

result= solver_DVM_3D(par);

%two different systems and so two rows
temp = cell(2,par.n_eqn);

for i = 1 : 2
    for j = 1 : par.n_eqn
        temp{i,j} = result(i,j).sol;
    end
end

density = compute_density(temp,par.Ax,par.Ay,par.all_w);
[ux,uy] = compute_velocity(temp,par.Ax,par.Ay,par.all_w);
theta = compute_theta(temp,par.Ax,par.Ay,par.all_w);

filename = 'result_Inflow2D_DVM.txt';
dlmwrite(filename,result(1,1).X','delimiter','\t','precision',10);
dlmwrite(filename,density','delimiter','\t','-append','precision',10);
dlmwrite(filename,ux','delimiter','\t','-append','precision',10);
dlmwrite(filename,uy','delimiter','\t','-append','precision',10);
dlmwrite(filename,theta','delimiter','\t','-append','precision',10);
end


% working on inflow boundaries, we consider vacum boundary conditions
function f = bc_inhomo(B,bc_id,Ax,Ay,id_sys,t)

    rho = 0;
    ux = 0;
    uy = 0;
    
    if t <= 1
        thetaIn = exp(-1/(1-(t-1)^2)) * exp(1);
    else
        thetaIn = 1;
    end
  
    
    f = diag(B) * 0;
    
%     id = find(diag(B) == 1);
%     
%     switch bc_id
%         % boundary at x = x_start
%         case 1
%             thetaIn = 0;
%     
%             for i = 1 : length(id)
%                 f(id(i)) = compute_fM(Ax,Ay,rho,ux,uy,thetaIn,id(i),id_sys);
%             end
%             
%         case 2
%             thetaIn = 0;
%     
%             for i = 1 : length(id)
%                 f(id(i)) = compute_fM(Ax,Ay,rho,ux,uy,thetaIn,id(i),id_sys);
%             end
%     end

end
        
% we have two systems of the same type 
function f = ic(x,id,Ax,Ay,sys_id)

density = exp(-(x-0.5).*(x-0.5)*100);
%density = zeros(length(x),1);

% velocity in the x-direction
vx = Ax(id,id);

% velocity in the y-direction
vy = Ay(id,id);

switch sys_id
    % value for g
    case 1  
        f = density * f0(vx,vy);
        
    % value for h
    case 2
        f = density * 0;

end

end

% f0 corresponding to the 2d velocity space
function f = f0(vx,vy)
f = exp(-(vx^2+vy^2)/2)/(2 * pi);
end

% compute the density on the complete grid
function f = compute_density(U,Ax,Ay,all_weights)

grid_points = length(U{1,1});

f = zeros(grid_points,1);

% loop over the grid
for i = 1 : grid_points
    
    % extract the value of g at a particular space point
    temp = cell2mat(cellfun(@(a) a(i),U(1,:),'Un',0));
    
    % might need to reshape
    if size(temp,1) ~= size(all_weights,1)
        temp = reshape(temp,size(all_weights));
    end
    
    f(i) = sum(all_weights.*temp);
end

end

% compute velocity in the x and y direction
function [vx,vy] = compute_velocity(U,Ax,Ay,all_weights)
grid_points = length(U{1,1});

vx = zeros(grid_points,1);
vy = zeros(grid_points,1);

% scale by velocity
all_weights_x = diag(Ax).*all_weights;
all_weights_y = diag(Ay).*all_weights;

% loop over the grid
for i = 1 : grid_points
    temp = cell2mat(cellfun(@(a) a(i),U(1,:),'Un',0));
    if size(temp,1) ~= size(all_weights_x,1)
        temp = reshape(temp,size(all_weights_x));
    end
    vx(i) = sum(all_weights_x.*temp);
    vy(i) = sum(all_weights_y.*temp);
end

end

% compute the temperature
function f = compute_theta(U,Ax,Ay,all_weights)
grid_points = length(U{1,1});

f = zeros(grid_points,1);

% scale by the He2(xi_1)
all_weights_x = (diag(Ax).*diag(Ax)-1).*all_weights/sqrt(2);

% scale by He2(xi_2)
all_weights_y = (diag(Ay).*diag(Ay)-1).*all_weights/sqrt(2);

% loop over the grid
for i = 1 : grid_points
    
    temp_g = cell2mat(cellfun(@(a) a(i),U(1,:),'Un',0));
    temp_h = cell2mat(cellfun(@(a) a(i),U(2,:),'Un',0));
    
    if size(temp_g,1) ~= size(all_weights_x,1)
        temp_g = reshape(temp_g,size(all_weights_x));
    end
    
    if size(temp_h,1) ~= size(all_weights_x,1)
        temp_h = reshape(temp_h,size(all_weights_x));
    end
    
    f(i) = sqrt(2) * sum((all_weights_x + all_weights_y).*temp_g + ...
                        all_weights.*temp_h)/3;
end

end

% compute the Maxwellian, there are two maxellians, one corresponding to
% g and one for h
function f = compute_fM(Ax,Ay,rho,ux,uy,theta,id,id_sys)

vx = Ax(id,id);
vy = Ay(id,id);


switch id_sys
    
    % maxwellian for g
    case 1
        f = (vx * ux + vy * uy + ((vx^2 + vy^2)/2 - 1) * theta + rho) * f0(vx,vy);
    % maxwellian for h
    case 2
        f = theta * f0(vx,vy)/sqrt(2);
end

end
