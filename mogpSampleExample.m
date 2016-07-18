%% GP sampling.
%% Clean
close all; clear all; clc;

%% Prepare kernel functions
Ksq = @(x,y,s,l) s*s*exp(-(x-y)^2/(2*l*l));
Ksin = @(x,y,s,l, p) s*s*exp(-2*sin(pi*(x-y)/p)^2/(l*l));

s1 = 1;
l1 = 0.1;
k1 = @(x,y) Ksq(x,y,s1,l1);

s2 = 0.5;
l2 = 10;
p2 = 2*pi;
k2 = @(x,y) Ksin(x,y,s2,l2, p2);


%% Sampling space
X = 0:0.01:5;
n = length(X);

%% Compute kernel
B1 = [1 0.5; 0.5 1.5];

K1 = [];
for i=1:n
    for j =1:n
        K1(i,j) = k1(X(i),X(j));
    end
end

B2 = [1 0.5; 0.5 1.5]
K2 = [];
for i=1:n
    for j =1:n
        K2(i,j) = k2(X(i),X(j));
    end
end


K = kron(B1,K1) + kron(B2,K2);

figure(1);
imagesc(K);
colorbar;


%% Sample
figure;
hold on;
u = randn(n*2,1);
% L  = chol(K); not very stable numerically 
[A S D] = svd(K);
L = A*sqrt(S);
F = L*u;
F = reshape(F,[n,2]);
plot(F(:,1));
plot(F(:,2));
