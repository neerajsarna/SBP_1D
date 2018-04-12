
function[Ax,B] = get_system_data(filenames)

%% boundary matrices for x = 1
B = create_sp_mat(filenames.B);


%% System matrix
Ax = create_sp_mat(filenames.Ax);
Ax = full(Ax);

end