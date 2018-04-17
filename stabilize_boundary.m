function BStable = stabilize_boundary(Ax,B)
n_eqn = size(Ax,2);

id_odd = 2:2:n_eqn;
id_even = 1:2:n_eqn;


% remove the highest order even moment
Trun_id_even = id_even(1:length(id_odd));

hatAoe = Ax(id_odd,Trun_id_even);

% Onsager matrix, ignored the minus
R = -B(:,Trun_id_even) * inv(hatAoe);

D = eig(R);

if ~isempty(find(D <0, 1))
    error('Onsager matrix not spd');
end
    
    
BStable = B;

% change all the even variables
BStable(:,id_even) = -R * Ax(id_odd,id_even);

end