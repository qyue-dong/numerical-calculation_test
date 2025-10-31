clear;
clc;
close all; 

[U1, u1, E1] = ErrorSolution(0.5, 0.25, 1/3,0.5);
[U2, u2, E2] = ErrorSolution(0.5, 0.25, 1/3,0.1);
[U3, u3, E3] = ErrorSolution(0.5, 0.25, 1/3,1);
x = -50:0.5:50-0.5;
U1=U1(:,21);
U2=U2(:,21);
U3=U3(:,21);
figure
plot(x,U1)
hold on
legend('c=0.5');
xlabel('x');
ylabel('$${u}$$(x,10)', 'Interpreter', 'LaTeX');
figure
plot(x, U2);
hold on
legend('c=0.1');
xlabel('x');
ylabel('$${u}$$(x,10)', 'Interpreter', 'LaTeX');
figure
plot(x, U3);
legend('c=1');
xlabel('x');
ylabel('$${u}$$(x,10)', 'Interpreter', 'LaTeX');