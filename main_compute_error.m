clear all;

M_values = 4:1:100;
grid_points = 300;
n_ref = 200;
t_end = 0.3;
filename = 'inflow';
foldername = 'result_Inflow';

[convg_rate,expected_rate] = compute_error(M_values,grid_points,n_ref,t_end,filename,foldername);

