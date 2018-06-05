
filename = '../result_Comp_DVM_Mom/result_HC2D_DVM_theta1_Kn0p1.txt';
result_DVM = dlmread(filename,'\t');

filename = '../result_HC2D/hc_tend_1_points_300_neqn_55.txt';
result_Mom = dlmread(filename,'\t');

theta = sqrt(2) * (result_Mom(4,:) + result_Mom(5,:) + result_Mom(6,:))/3;

figure(1);
plot(result_DVM(1,:),result_DVM(5,:),'-r',result_Mom(1,:),theta,'-k');
grid on;
legend('dvm','mom');
