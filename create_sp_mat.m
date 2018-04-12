% given the filename create the sp matrix
function[sp_mat] = create_sp_mat(filename)

    fileID = fopen(filename,"r");
    indices = dlmread(filename);

    ia = indices(2:end,1);
    ja = indices(2:end,2);
    v = indices(2:end,3);

    % CAUTION: THE NUMBER OF NON-ZEROS NOT READ, RATHER THE DIMENSIONS ARE
    % BEING READ.
    m = indices(1,1);
    n = indices(1,2);

    sp_mat = sparse(ia,ja,v,m,n);

    fclose('all');
end