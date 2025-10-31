clear;
clc;
close all; 

[U, u, E1] = ErrorSolution(0.25, 0.25,1/3,1);
x = -50:0.25:50-0.25;
t = 0:0.25:10;
[x, t] = meshgrid(x, t);

% U(x,t) 三维表面图
figure('Color', 'w');
U=U';
mesh(x, t, U);
colormap(cool); % 你也可以试试 parula, hot, cool 等
colorbar;
xlabel('x', 'FontSize', 14);
ylabel('t', 'FontSize', 14);
zlabel('U(x,t)', 'FontSize', 14);
title('Numerical Solution', 'FontSize', 16);
set(gca, 'FontSize', 14);
view(45, 30); % 调整视角
shading interp; % 平滑色彩
figure('Color', 'w');
u=u';
mesh(x, t, u);
colormap(hsv); % 你也可以试试 parula, hot, cool 等
colorbar;
xlabel('x', 'FontSize', 14);
ylabel('t', 'FontSize', 14);
zlabel('u(x,t)', 'FontSize', 14);
title('Exact Solution', 'FontSize', 16);
set(gca, 'FontSize', 14);
view(45, 30); % 调整视角
shading interp; % 平滑色彩
error=abs(u-U);
figure('Color', 'w');
mesh(x, t, error);
colormap(hsv); % 你也可以试试 parula, hot, cool 等
colorbar;
xlabel('x', 'FontSize', 14);
ylabel('t', 'FontSize', 14);
zlabel('error(x,t)', 'FontSize', 14);
title('Error graph', 'FontSize', 16);
set(gca, 'FontSize', 14);
view(45, 30); % 调整视角
shading interp; % 平滑色彩