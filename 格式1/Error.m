clear;
clc;
close all;
[U,u,~]=ErrorSolution(0.5,0.1,1/3,0.1);
xl = -50;
xr = 50;
x=xl+0.5:0.5:xr;
t=0:0.25:10;
figure
plot(x, U(:,1), 'm-')
hold on
plot(x, U(:,11), 'b-')
hold on
plot(x, U(:,21), 'c-')
xlabel('x')
ylabel('u')
hold on
plot(x, U(:,31), 'r-')
hold on
plot(x, U(:,41), 'g-')
xlabel('x')
ylabel('u(x,t)')
legend({'t=0', 't=2.5', 't=5', 't=7.5', 't=10'}, ...
       'Location', 'northeast', ...  % 图例位置（右上角）
       'FontSize', 10);


