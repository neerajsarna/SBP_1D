clear all;

sol_dvm = dlmread('../result_Inflow/DVM_inflow_tend_0.3_points_300_neqn_50.txt');
sol_dvm_f = dlmread('../result_Inflow/DVM_inflow_0.3_points_50_neqn_20.txt');
sol_mom = dlmread('../result_Inflow/inflow_tend_0.3_points_300_neqn_200.txt');

initial_rho = exp(-(sol_mom(1,:)-0.5).*(sol_mom(1,:)-0.5)*100);

% plot for density
figure(1);
plot(sol_dvm(1,:),sol_dvm(2,:),'-*',sol_mom(1,:),sol_mom(2,:),'-o',sol_mom(1,:), ...
            initial_rho,'-o','markersize',4);
grid on;
legend('DVM','M=200','initial condition','Location','best');
title('variation of density');
ylabel('\mu_0(f)','FontSize',18);
xlabel('x','FontSize',18);
xt = get(gca, 'YTick');
set(gca, 'FontSize', 16);

% plot for velocity
figure(2);
plot(sol_dvm(1,:),sol_dvm(3,:),'-*',sol_mom(1,:),sol_mom(3,:),'-o','markersize',4);
grid on;
legend('DVM','M=200','Location','best');
title('variation of velocity');
ylabel('\mu_1(f)','FontSize',18);
xlabel('x','FontSize',18);
xt = get(gca, 'YTick');
set(gca, 'FontSize', 16);

% plot for temperature
figure(3);
plot(sol_dvm(1,:),sol_dvm(4,:) * sqrt(2),'-*',sol_mom(1,:),sol_mom(4,:) * sqrt(2),'-o','markersize',4);
grid on;
legend('DVM','M=200','Location','best');
title('variation of temperature');
h = ylabel('$\sqrt{2}\mu_2(f)$','FontSize',18);
set(h,'Interpreter','latex') 
xlabel('x','FontSize',18);
xt = get(gca, 'YTick');
set(gca, 'FontSize', 16);



[x_grid,v_grid] = meshgrid(sol_dvm_f(1,:),sol_dvm_f(2,1:size(sol_dvm_f,1)-2));

figure(4);
surf(x_grid,v_grid,sol_dvm_f(3:end,:));
xlabel('x','FontSize',18);
h = ylabel('$\xi$','FontSize',18);
set(h,'Interpreter','latex');
zlabel('f');
title('variation of kinectic solution at T = 0.3');
xt = get(gca, 'YTick');
set(gca, 'FontSize', 16);