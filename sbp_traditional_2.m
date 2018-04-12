function[D,P,InvP] = sbp_traditional_2(delta_x,N)
% N is the number of intervals so we have N + 1 grid points

% 2 for the boundary and we have two boundaries. Total N+1-2 = N-1 interior
% points
nz = 2 * 2 + 2 * (N-1);
IA = zeros(nz,1);
JA = zeros(nz,1);
VA = zeros(nz,1);

% left boundary
IA(1:2) = [1 1];
JA(1:2) = [1 2];
VA(1:2) = [-1 1];

% right boundary
IA(end-1 : end) = [N+1 N+1];
JA(end-1 : end) = [N N+1];
VA(end-1 : end) = [-1 1];

% loop over interior points
for i=1:N-1
    IA(2 * i + 1 : 2 * (i + 1)) = [i+1 i+1];
    JA(2 * i + 1 : 2 * (i + 1)) = [i i+2];
    VA(2 * i + 1 : 2 * (i + 1)) = [-1/2 1/2];
end

D = sparse(IA,JA,VA,N+1,N+1)/delta_x;

IA = zeros(N+1,1);
JA = zeros(N+1,1);
VA = zeros(N+1,1);

for i = 1 : N+1
    IA(i) = i;
    JA(i) = IA(i);
    VA(i) = 1;
end

VA(1) = 1/2;
VA(end) = 1/2;

P = delta_x * sparse(IA,JA,VA,N+1,N+1);

VA(1) = 2;
VA(end) = 2;

InvP = sparse(IA,JA,VA,N+1,N+1)/delta_x;


end

