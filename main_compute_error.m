
clear all;

%% routines for the 1D case
M_values = [4:2:24,5:2:25];
M_values = sort([M_values,5]);
grid_points = 300;
n_ref = 200;
n_ref2 = 200;
t_end = 0.3;
filename = 'wall';
foldername = 'result_collision_gaussian_1x1v';

compute_error(M_values,grid_points,n_ref,t_end,filename,foldername);

%% Heat Conduction

% M_values = sort([4:4:24,3:4:27]);
% grid_points = 300;
% n_ref = 55;
% t_end = 1;
% filename = 'hc';
% foldername = 'result_HC2D';
% 
% [convg_rate_Odd,convg_rate_Even,expected_rate] = compute_error2D(M_values,grid_points,n_ref, ...
%                                             t_end,filename,foldername);