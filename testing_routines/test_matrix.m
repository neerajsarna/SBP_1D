function test_matrix(n_eqn)

B1 = dvlp_BInflow1D(n_eqn);
A1 = dvlp_Ax1D(n_eqn);

B1 = stabilize_boundary(A1,B1);

end