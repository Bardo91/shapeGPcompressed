%% GP sampling.

%% Prepare kernel functions
Ksq = @(x,y,s,l) s*s*exp(-(x-y)^2/(2*l*l));
Ksin = @(x,y,s,l, p) s*s*exp(-2*sin(pi*(x-y)/p)^2/(l*l));

s1 = 0.2;
l1 = 0.5;
k1 = @(x,y) Ksq(x,y,s1,l1);

s2 = 0.5;
l2 = 1;
p2 = 1;
k2 = @(x,y) Ksin(x,y,s2,l2, p2);


%% Sampling space and data points
X = 0:0.01:5;
n = length(X);

fs = [0, 0.3, 0.1, -0.1 -0.2, 0]';
Xs = [0.4, 2.3,4.2];
m = length(Xs);

%% Compute kernel
B1 = [1.2 0.5; 0.5 1.5];

K1 = [];
for i=1:n
    for j =1:n
        K1(i,j) = k1(X(i),X(j));
    end
end
K1 = kron(B1,K1);


Kss1 = [];
for i=1:m
    for j =1:m
        Kss1(i,j) = k1(Xs(i),Xs(j));
    end
end
Kss1 = kron(B1,Kss1);

Ks1 = [];
for i=1:m
    for j =1:n
        Ks1(i,j) = k1(Xs(i),X(j));
    end
end
Ks1 = kron(B1,Ks1);

% %
% B2 = [1 0.5; 0.5 1];
% 
% K2 = [];
% for i=1:n
%     for j =1:n
%         K2(i,j) = k2(X(i),X(j));
%     end
% end
% K2 = kron(B2,K2);
% 
% 
% Kss2 = [];
% for i=1:m
%     for j =1:m
%         Kss2(i,j) = k2(Xs(i),Xs(j));
%     end
% end
% Kss2 = kron(B2,Kss2);
% 
% Ks2 = [];
% for i=1:m
%     for j =1:n
%         Ks2(i,j) = k2(Xs(i),X(j));
%     end
% end
% Ks2 = kron(B2,Ks2);


K = K1;%+K2;
Ks = Ks1;%+Ks2;
Kss = Kss1;%+Kss2;

%% Regression
nus = 0 + Ks'*(Kss\fs);
sig = K - Ks'*(Kss\Ks);

nus = reshape(nus,[n,2]);
d = diag(sig);
d = reshape(d,[n,2]);
fs = reshape(fs,[m,2]);

figure;
hold on;
plot(X, nus(:,1),'r');
plot(Xs,fs(:,1), 'r*');
plot(X, nus(:,1)+real(sqrt(d(:,1))),'k');
plot(X, nus(:,1)-real(sqrt(d(:,1))),'k');


plot(X, nus(:,2),'b');
plot(Xs,fs(:,2), 'b*');
plot(X, nus(:,2)+real(sqrt(d(:,2))),'k');
plot(X, nus(:,2)-real(sqrt(d(:,2))),'k');