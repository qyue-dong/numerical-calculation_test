
clear;
clc;
close all; 
[U, u,E1]=ErrorSolution(0.1, 0.01, 1/3,1); 
U1=U(:,end)-u(:,end);
error1=sqrt(2* sum(U1.^2));
error2=norm(U1,inf);

figure
t=0:0.01:10;
plot(t, E1);
x1=E1(:,1)
x2=E1(:,201)
x3=E1(:,401)
x4=E1(:,601)
x5=E1(:,801)
x6=E1(:,1001)
%ylim([5.599999,5.6])
xlabel('t');
ylabel('$$\tilde{E}_1^n$$', 'Interpreter', 'LaTeX');
