close all; clear; clc
%% generate synthetic data
N = 1000; T = 1000;
p = 4 ; rho = 0.5;
theta = 100; snr = 1/theta;% control the snr
e = zeros(N,T);
for i=1:N
    v = randn(1,T);
    e(i,1) = randn();
    for j = 2:T
        e(i,j) = rho*e(i,j-1)+v(j);
    end
end
u = sqrt(1-rho^2)*e;
L = randn(N,p);
F = randn(p,T);
R = L*F+sqrt(theta*p)*u;
%% real eigenvalue of the covariance matrix of p-level residual
p_real = p;
E = eig_real(R,p_real);  % compute the real eigenvalue
%% comparison between the ESD and ASD
a = load('data');
plot(0.01:0.1:9.01,a.pr,'r','LineWidth',1);
hold on
Nbins = 50;
[Y,X] = hist(E,Nbins);
r = max(E)-min(E);
Y = Y/sum(Y)*Nbins/r;
bar(X,Y,'FaceColor','none')
pic_title = strcat('snr=',num2str(snr),',numbers of removed factors = ',num2str(p_real));
% pic_title = strcat('snr=0',',numbers of removed factors = ',num2str(p_real));
title(pic_title)

%% RMSE
u = mapstd(u);
c1=u*u'/T;
E0 = eig(c1);
RMSE = sum((E -E0).^2)/N; 