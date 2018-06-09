clear all;

sol_dvm = dlmread('../result_collision_gaussian_1x1v/DVM_wall_tend_0.3_points_50_neqn_50.txt');
sol_mom = dlmread('../result_collision_gaussian_1x1v/wall_tend_0.3_points_50_neqn_5.txt');
sol_mom2 = dlmread('../result_collision_gaussian_1x1v/wall_tend_0.3_points_50_neqn_25.txt');
sol_mom3 = dlmread('../result_collision_gaussian_1x1v/wall_tend_0.3_points_50_neqn_100.txt');

% plot for density
figure(1);
plot(sol_dvm(1,:),sol_dvm(2,:),'-*',...
     sol_mom(1,:),sol_mom(2,:),'-o', ...
     sol_mom2(1,:),sol_mom2(2,:),'-o',...
      sol_mom3(1,:),sol_mom3(2,:),'-o','markersize',4);
  
grid on;
legend('DVM','M=5','M=25','M=70','Location','best');
title('variation of density');
ylabel('\mu_0(f)','FontSize',18);
xlabel('x','FontSize',18);
xt = get(gca, 'YTick');
set(gca, 'FontSize', 16);

% plot for velocity
figure(2);
plot(sol_dvm(1,:),sol_dvm(3,:),'-*',...
    sol_mom(1,:),sol_mom(3,:),'-o',...
    sol_mom2(1,:),sol_mom2(3,:),'-o',...
    sol_mom3(1,:),sol_mom3(3,:),'-o',...
    'markersize',4);
grid on;
legend('DVM','M=5','M=25','M=70','Location','best');
title('variation of velocity');
ylabel('\mu_1(f)','FontSize',18);
xlabel('x','FontSize',18);
xt = get(gca, 'YTick');
set(gca, 'FontSize', 16);

% plot for temperature
figure(3);
plot(sol_dvm(1,:),sol_dvm(4,:) * sqrt(2),'-*',...
    sol_mom(1,:),sol_mom(4,:) * sqrt(2),'-o',...
    sol_mom2(1,:),sol_mom2(4,:) * sqrt(2),'-o',...
    sol_mom3(1,:),sol_mom3(4,:) * sqrt(2),'-o',...
    'markersize',4);
grid on;
legend('DVM','M=5','M=25','M=70','Location','best');
title('variation of temperature');
h = ylabel('$\sqrt{2}\mu_2(f)$','FontSize',18);
set(h,'Interpreter','latex') 
xlabel('x','FontSize',18);
xt = get(gca, 'YTick');
set(gca, 'FontSize', 16);

%% comparison of distribution function
% sol_dvm_f = dlmread('../result_Inflow_fluctuateT_1x1v/DVM_f_inflow_1_points_300_neqn_50.txt');
% [x_mesh,v_mesh] = meshgrid(sol_dvm_f(1,:),sol_dvm_f(2,1:size(sol_dvm_f,1)-2));
% 
% figure(4)
% surf(x_mesh,v_mesh,sol_dvm_f(3:end,:));
% title('variation of f');
