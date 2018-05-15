
clear all;

%% routines for the 1D case
% M_values = 4:1:40;
% grid_points = 300;
% n_ref = 200;
% n_ref2 = 199;
% t_end = 0.3;
% filename = 'inflow';
% foldername = 'result_Inflow';
% 
% [convg_rate,expected_rate] = compute_error(M_values,grid_points,n_ref,n_ref2,t_end,filename,foldername);

%% Heat Conduction

M_values = sort([4:4:24,3:4:27]);
grid_points = 300;
n_ref = 55;
t_end = 1;
filename = 'hc';
foldername = 'result_HC2D';

[convg_rate_Odd,convg_rate_Even,expected_rate] = compute_error2D(M_values,grid_points,n_ref, ...
                                            t_end,filename,foldername);