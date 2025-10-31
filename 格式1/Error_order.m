clear;
clc;
close all;

[U1,u1,~]=ErrorSolution3(0.5,0.25,1/3,1);
[U2,u2,~]=ErrorSolution3(0.25,0.0625,1/3,1);
[U3,u3,~]=ErrorSolution3(0.125,0.015625,1/3,1);
[U4,u4,~]=ErrorSolution3(0.0625,0.00390625,1/3,1);
U2=U2(2:2:end,1:4:end);
u2=u2(2:2:end,1:4:end);
U3=U3(4:4:end,1:16:end);
u3=u3(4:4:end,1:16:end);
U4=U4(8:8:end,1:64:end);
u4=u4(8:8:end,1:64:end);
error1=u1(:,end)-U1(:,end);
error2=u2(:,end)-U2(:,end);
error3=u3(:,end)-U3(:,end);
error4=u4(:,end)-U4(:,end);
%...........无穷范数
%tif1=norm(error1,inf);
%tif2=norm(error2,inf);
%tif3=norm(error3,inf);
%tif4=norm(error4,inf);
%。。。。。。二范数
tif1=sqrt(0.5 * sum(error1.^2));
tif2=sqrt(0.5 * sum(error2.^2));
tif3=sqrt(0.5 * sum(error3.^2));
tif4=sqrt(0.5 * sum(error4.^2));
result=log2(tif1/tif2);
result1=log2(tif2/tif3);
result2=log2(tif3/tif4);





