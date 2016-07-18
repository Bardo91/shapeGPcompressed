%% Clean workspace
clear all; close all; clc;

%% Generate Data
A = rand(2);
m = 10;
angles = linspace(0,2*pi,m);
f = @(theta) [cos(theta),sin(theta)];
r = arrayfun(@(i)f(angles(i)),1:m,'UniformOutput',false);
X = reshape(cell2mat(r),[2 m]);
norms = X;
angles = tan(norms(2,:)./norms(1,:));
m = length(X);


f = zeros(m,1);
f = [f,angles'];
[a b] = size(f);
f = reshape(f', a*b, 1);

%% Kernels
sigma1 = 1;
length1 = 0.5;
k1 = @(x,y)sigma1^2 * exp(norm(x-y)^2/(2*length1^2));

sigma2 = 1;
length2 = 0.5;
p = 2*pi;
k2 = @(a,b)sigma2^2 * exp(2*sin(pi*norm(a-b)/p)^2/(lengh2^2));

%%  Create Evaluation space
limx = 2;
limy = 2;
step = 0.1;
[Xg,Yg] = meshgrid(-limx:step:limx,-limy:step:limy);
[d1,d2] = size(Xg);
Xs = [reshape(Xg,d1*d2,1),reshape(Yg,d1*d2,1)]';
n = length(Xs);

%% Evaluate kernels
K = zeros(2*m,2*m);
for i=1:m
    for j=1:m
        K(i,j)=k1();
        K(i,j+1)=;
        K(i+1,j)=;
        K(i+1,j+1)=;
    end
end

Ks = zeros(2*n,2*m);
for i=1:n
    for j=1:m
        Ks(i,j)=;
        Ks(i,j+1)=;
        Ks(i+1,j)=;
        Ks(i+1,j+1)=;
    end
end

Kss = zeros(2*n,2*n);
for i=1:n
    for j=1:n
        Kss(i,j)=;
        Kss(i,j+1)=;
        Kss(i+1,j)=;
        Kss(i+1,j+1)=;
    end
end


%% Computing means
% display('Computing means');
% mu = zeros(m*4,1);
% for i = 1:m
%     mu((i-1)*4 +1) = meanValue(X(:,i));
%     mu((i-1)*4 +2:(i-1)*4 +4) = meanGrad(X(:,i));
% end
% mus = zeros(n*4,1);
% for i = 1:n
%     mus((i-1)*4 +1) = meanValue(Xs(:,i));
%     mus((i-1)*4 +2:(i-1)*4 +4) = meanGrad(Xs(:,i));
% end

%% Regression
display('Computing regression');
R = K\Ks;
fs = mus + R'*(f - mu);
var = diag(Kss) - diag(R'*Ks);

sig = sqrt(var);
SIG = reshape(sig(1:4:end), d1,d2,d3);
Fs = reshape(fs(1:4:end),d1,d2,d3);

figure
contourf(Xg,Yg,Fs);        % draw image and scale colormap to values range
colorbar;
axis equal
