close all;clear;clc
T = 2000;
N = 1000;
p = 4;
%% generate synthetic data
L = randn(N,p);
F =randn(p,T);
U = randn(N,T); %white noise
snr = 95;
R =snr* L*F+U;
p_real = 4;
E0 = eig_real(R,p_real);

[~,T] = size(R);
[F,~] = eigs(R'*R,p_real); % obtain the p largest eigenvalues and corresponding eigenvectors
F=F';
L =  R*F'/(F*F');
U1 = R-L*F;

%%
Nbins = 100;
[Y,X] = hist(E0,Nbins);
r = max(E0)-min(E0);
Y = Y/sum(Y)*Nbins/r;
bar(X,Y,'FaceColor','none')
hold on
%%
c1 = N/T;
a = (1-sqrt(c1))^2; b=(1+sqrt(c1))^2;
t = a:0.01: b+0.01;
u = ((b-t).*(t- a)).^(1/2)./(2*pi*c1.*t);
plot(t,u,'r','LineWidth',1);
pic_title = strcat('numbers of removed factors = ',num2str(p_real));
title(pic_title)
% save 
