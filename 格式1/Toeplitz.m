function [T] = Toeplitz(n, first, second, third)
t=zeros(1,n);
size=length(t);
t(1)=first;
t(2)=second;
t(size)=third;
New=ones(1,size);
New(1)=t(1);
New(2:size)=t(size:-1:2);
T=toeplitz(New,t)
end
