% compute the eigenvector corresponding to the negative eigenvalues
function Xminus = compute_Xminus(A)
[V,D] = eig(A);
D= D(sub2ind(size(D),1:size(D,1),1:size(D,2)));
loc_neg = find(D<-(1e-10));
Xminus = V(:,loc_neg);
end
