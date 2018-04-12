clear all;

%========================================================================
% Problem Parameters
%========================================================================
par = struct(...
'name','Lattice Test',... % name of example
'n_mom',3,... % order of moment approximation
'sigma_a',@sigma_a,... % absorption coefficient (defined below)
'sigma_s0',@sigma_s0,... % isotropic scattering coefficient (def. below)
'source',@source,... % source term (defined below)
'ax',[0 7 0 7],... % coordinates of computational domain
'n',[100 100],... % numbers of grid cells in each coordinate direction
 't_end',1.0,... % the end time of the computatio
 'diff_order',4,... % the difference order in the physical space
 'CFL',2.0,...      % the crude cfl number
 'num_bc',4,... % number of boundaries in the domain
'output',@output... % problem-specific output routine (defined below
);

par.t_plot = linspace(0,par.t_end,50);
   
par.n_eqn = (par.n_mom + 1) * (par.n_mom + 2)/2;

[par.system_data.Ax,par.system_data.Ay,par.mom_order] = closure_pn(par.n_mom);

% the moment variable to be output
par.mom_output = 1;

% create the penalty matrix and penalty * B for all the boundaries
% a loop over all the boundaries
% for i = 1 : par.num_bc
%     alpha = (i-1) * pi/2;
%     projector =  global_projector(par.n_mom,alpha, ...
%                                   par.system_data.Perm,par.system_data.InvPerm);
%                               
%     par.system_data.penalty_B{i}  =  full(projector' * par.system_data.Sigma * ...
%                                      par.system_data.BInflow * projector);
%                                  
%     par.system_data.penalty{i}  =  full(projector' * par.system_data.Sigma);
% end

% solve the system
solution = solver(par);

%========================================================================
% Problem Specific Functions
%========================================================================
function f = sigma_a(x,y)
% Absorption coefficient.
cx = ceil(x); cy = ceil(y);
g = (ceil((cx+cy)/2)*2==(cx+cy)).*(1<cx&cx<7&1<cy&cy-2*abs(cx-4)<4);
f = (1-g)*0+g*10;
end

function f = sigma_s0(x,y)
% Isotropic scattering coefficient.
cx = ceil(x); cy = ceil(y);
g = (ceil((cx+cy)/2)*2==(cx+cy)).*(1<cx&cx<7&1<cy&cy-2*abs(cx-4)<4);
f = (1-g)*1+g*0;
end

function f = source(x,y)
% Radiation source (only for zeroth moment).
f = 3<x&x<4&3<y&y<4;
end

function output(par,x,y,U,step)
% Plotting routine.
cax = [-7 0];                    % Colormap range used for log10 scaling.
vneg = cax(1)-diff(cax)/254;        % Value assigned where U is negative.
V = log10(max(U,1e-50));% Cap U s.t. U>0 and use logarithmic color scale.
Vcm = max(V,cax(1));                      % Cap colormap plot from below.
Vcm(U<0) = vneg;              % Assign special value where U is negative.
clf, subplot(1,3,1:2)
imagesc(x,y,Vcm'), axis xy equal tight, caxis([vneg cax(2)])
title(sprintf('t = %0.2f',par.t_plot(step)))
cm = jet(256); cm(1,:) = [1 1 1]*.5; % Change lowest color entry to gray.
colormap(cm), colorbar('ylim',cax)
xlabel('x'), ylabel('y')
subplot(1,3,3)
plot(y,interp2(x,y,V',y*0+3.5,y))   % Evaluate solution along line x=3.5.
axis([par.ax(1:2) cax+[-1 1]*.1])
title('Cut at x=3.5'), xlabel('y')
drawnow
end


function[filenames] = dvlp_filenames(nEqn)

filenames = struct;

filenames.BInflow = strcat("generic_2D/system_mat/BInflow_",num2str(nEqn));
filenames.BInflow = strcat(filenames.BInflow,".txt");

filenames.Perm = strcat("generic_2D/system_mat/Perm_",num2str(nEqn));
filenames.Perm = strcat(filenames.Perm,".txt");

filenames.InvPerm = strcat("generic_2D/system_mat/InvPerm_",num2str(nEqn));
filenames.InvPerm = strcat(filenames.InvPerm,".txt");

filenames.Ax = strcat("generic_2D/system_mat/Ax_",num2str(nEqn));
filenames.Ax = strcat(filenames.Ax,".txt");

filenames.Ay = strcat("generic_2D/system_mat/Ay_",num2str(nEqn));
filenames.Ay = strcat(filenames.Ay,".txt");

filenames.Sigma = strcat("generic_2D/system_mat/Sigma_",num2str(nEqn));
filenames.Sigma = strcat(filenames.Sigma,".txt");
end
