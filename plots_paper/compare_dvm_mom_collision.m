clear all;

% compatible and incompatible boundary conditions
sol_dvm_f = dlmread('../result_collision_gaussian_1x1v/DVM_f_wall_tend_0.3_points_300_neqn_100.txt');
sol_dvm = dlmread('../result_collision_gaussian_1x1v/DVM_wall_tend_0.3_points_300_neqn_100.txt');
sol_mom = dlmread('../result_collision_gaussian_1x1v/wall_tend_0.3_points_300_neqn_200.txt');


valuesHermite = dlmread('../result_collision_gaussian_1x1v/valuesHermite.txt');
xPoints = sol_dvm_f(1,:);
vPoints = sol_dvm_f(2,1:(length(sol_dvm_f(:,1))-3));

f_mom = valuesHermite * sol_mom(2:end,:);

f_mom_Even = valuesHermite(:,1:2:end) * sol_mom(2:2:end,:);
f_mom_Odd = valuesHermite(:,2:2:end) * sol_mom(3:2:end,:);

error = log10(abs(f_mom - sol_dvm_f(4:end,:)));

% meshlines = 10;
% plot_surf(xPoints,vPoints,sol_dvm_f(4:end,:),meshlines,1);
% plot_surf(xPoints,vPoints,f_mom,meshlines,2);
% plot_surf(xPoints,vPoints,error,meshlines,3);
% 
% figure(4);
% plot(sol_mom(1,:),sol_mom(2,:),'-o', ...
%     sol_dvm(1,:),sol_dvm(2,:),'-*', ...
%      'markersize',4);
% grid on;
% legend('M=200','DVM','Location','best');
% title('variation of deviation in density');
% xlabel('x','FontSize',18);
% xt = get(gca, 'YTick');
% set(gca, 'FontSize', 16);


%% comparing the distribution function at the wall
figure(5);
plot(vPoints,f_mom(:,1),'-o', ...
    vPoints,sol_dvm_f(4:end,1),'-*', ...
     'markersize',4);
grid on;
legend('M=200','DVM','Location','best');
h = xlabel('$\xi$','FontSize',18);
set(h,'Interpreter','latex');
title('kinetic solution x = 0');
xt = get(gca, 'YTick');
xlim([-4,4]);
set(gca, 'FontSize', 16);

figure(6);
plot(vPoints,f_mom(:,10),'-o', ...
    vPoints,sol_dvm_f(4:end,10),'-*', ...
     'markersize',4);
grid on;
legend('M=200','DVM','Location','best');
h = xlabel('$\xi$','FontSize',18);
set(h,'Interpreter','latex');
title('kinetic solution x = 1/30');
xlim([-4,4]);
set(gca, 'FontSize', 16);

figure(7);
plot(vPoints,f_mom_Odd(:,1),'-o', ...
    vPoints,f_mom_Even(:,1),'-*', ...
     'markersize',4);
grid on;
legend('Odd Part(M=200)','Even Part(M=200)','Location','best');
h = xlabel('$\xi$','FontSize',18);
set(h,'Interpreter','latex');
title('kinetic solution at x = 0');
xt = get(gca, 'YTick');
xlim([-4,4]);
set(gca, 'FontSize', 16);

%% comparing the density profiles 
sol_mom_M1 = dlmread('../result_collision_gaussian_1x1v/wall_tend_0.3_points_300_neqn_25.txt');
sol_mom_M2 = dlmread('../result_collision_gaussian_1x1v/wall_tend_0.3_points_300_neqn_10.txt');
sol_mom_M3 = dlmread('../result_collision_gaussian_1x1v/wall_tend_0.3_points_300_neqn_5.txt');

% plot for density
figure(1);
plot(sol_mom(1,:),sol_mom(2,:),'-o', ...
    sol_mom_M1(1,:),sol_mom_M1(2,:),'-*', ...
     sol_mom_M2(1,:),sol_mom_M2(2,:),'-^', ...
     sol_mom_M3(1,:),sol_mom_M3(2,:),'-<', ...
     'markersize',4);
grid on;
legend('M=200','M=25','M=10','M=5','Location','best');
title('variation of deviation in density');
xlabel('x','FontSize',18);
xt = get(gca, 'YTick');
set(gca, 'FontSize', 16);

function f = plot_surf(x,y,z,meshlines,plot_id)
%% Changing(reducing) the number of lines (edges / mesh) shown by a surface plot.
figure(plot_id);
s = surf(x,y,z,...
     'FaceColor','interp','EdgeColor','none');
xlabel('x','FontSize',18)
h = ylabel('$\xi$','FontSize',18);
set(h,'Interpreter','latex');
zlabel('f')
xt = get(gca, 'YTick');
set(gca, 'FontSize', 16);

%% Extract X,Y and Z data from surface plot
x=s.XData;
y=s.YData;
z=s.ZData;

% For R2014a and earlier:
% x=get(s,'XData');
% y=get(s,'YData');
% z=get(s,'ZData');


%% Divide the lengths by the number of lines needed
xnumlines = meshlines; % 10 lines
ynumlines = meshlines; % 10 partitions
xspacing = round(length(x)/xnumlines);
yspacing = round(length(y)/ynumlines);

%% Plot the mesh lines 
% Plotting lines in the X-Z plane
hold on

for i = 1:yspacing:length(y)
    Y1 = y(i)*ones(size(x)); % a constant vector
    Z1 = z(i,:);
    plot3(x,Y1,Z1,'-k');
end

% Plotting lines in the Y-Z plane

for i = 1:xspacing:length(x)
    X2 = x(i)*ones(size(y)); % a constant vector
    Z2 = z(:,i);
    plot3(X2,y,Z2,'-k');
end

hold off


end

