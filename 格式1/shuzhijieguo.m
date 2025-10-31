clear;
clc;
close all;
[U,u,~]=ErrorSolution(0.5,0.25,1/3,1);
xl = -50;
xr = 50;
x=xl+0.5:0.5:xr;
t=0:0.25:10;
figure
plot(x, U(:,11), 'r-.^')
legend('t=5')
hold on
plot(x, U(:,21), 'b-.^')
legend('t=10')
hold on
plot(x, U(:,31), 'c-.^')
legend('t=15')
hold on
plot(x, U(:,41), 'g-.^')
legend('t=2.5', 't=5', 't=7.5', 't=10')
xlabel('x')
ylabel('u')