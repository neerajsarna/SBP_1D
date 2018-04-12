function modA = compute_Amod(A)
[V,D] = eig(A);
modA = V * abs(D)/V;
end

