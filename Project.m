clear all
close all
clc

N = 101;
u = randi([0 1], N, 1);
ind = find(u == 0);
u(ind) = -1;
mu  = 0.97;
T = zeros(N+1,N+1);
T(1,1) = sqrt(mu)*u(1);
for i = 2 : N
    T(i,1) = sqrt(mu)*u(i);
    v      = zeros(1,N+1);
    for k = 1 : i-1
        s = zeros(1,N+1);
        s(1,k+1) = 1;
        v  = v - mu*u(k)*u(i).*[T(k,:) + s];
    end
    T(i,:) = T(i,:) + v;
end
%% 

mu  = 0.97;
T1 = zeros(N+1,N+1);
T1(1,1) = sqrt(mu)*u(1);
p = zeros(1,N+1);
p(1) = mu;
for i = 2 : N
    T1(i,1) = sqrt(p(1)/(1+p(1)))*u(i);
    v      = zeros(1,N+1);
    p(i) = p(i-1)/(1+p(i-1));
    for k = 1 : i-1
        s = zeros(1,N+1);
        s(1,k+1) = 1;
        v  = v - p(i-1)/(1+p(i-1))*u(k)*u(i).*[T1(k,:) + s];
    end
    T1(i,:) = T1(i,:) + v;
end

