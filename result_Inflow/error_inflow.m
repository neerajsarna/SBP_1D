clear all;

% we compute the error from 
% number of intervals
n = 150;

% spacing between two grid points
delta_x = 1/150;

% end time
t_end = 0.3;

% we need the P matrix for computing the error
[D,P,InvP] = sbp_traditional_2(delta_x,n);

num_samples = 27;

result = cell(num_samples,1);
counter = 1;

% loop over all the systems 
for i = 4:30
    output_filename = strcat('inflow_tend_',num2str(t_end),'_points_',num2str(n),'_neqn_');
    output_filename = strcat(output_filename,num2str(i),'.txt');
    result{counter} = dlmread(output_filename,'\t');
    counter = counter + 1;
end

errorTot = [];

% computing the error 
for i = 1 : num_samples-1
    
    % num rows we need to extract from the higher order moment solution
    num_rows = size(result{i},1);
    
    error = [(result{num_samples}(1:num_rows,:)-result{i}); result{num_samples}(num_rows+1:end,:)];
    
    % we need to transpose to allow multiplication with P
    error = error';

    % compute the error. The quantity dot(error,P*error,1) gives us ||.||^2_{L^2(\Omega;\mbb R^n)}
    errorTot = [errorTot sqrt(sum(dot(error,P*error,1)))];
    
end

loglog(errorTot);
