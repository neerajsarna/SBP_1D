
clear all;

%% routines for the 1D case
% M_values = [4:2:40,44:2:70,5:2:41,45:2:71];
% M_values = [44:4:70,45:4:71,75:5:105];
% M_values = 45:5:105;
% M_values = sort(M_values);
% grid_points = 300;
% n_ref = 200;
% n_ref2 = 200;
% t_end = 0.3;
% filename = 'wall';
% foldername = 'result_collision_gaussian_1x1v';
% 
% error = compute_error(M_values,grid_points,n_ref,t_end,filename,foldername);

%% pulsating temp

M_values = sort([4:4:24,3:4:27]);
grid_points = 300;
n_ref = 45;
n_ref2 = 44;
t_end = 0.5;
filename = 'inflow';
foldername = 'result_Inflow_1x3v_Kn0p1_cos_theta1';

[convg_rate_Odd,convg_rate_Even,expected_rate] = compute_error2D(M_values,grid_points,n_ref,n_ref2, ...
                                            t_end,filename,foldername);