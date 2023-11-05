% test_cap_spectral
% 29-Aug-2005
% Q.T. Le Gia, UNSW, Sydney, Australia.
clear all;

% we want solution on a spherical cap centered at ctradius alpha
Nr = 10000;
ct = [1/sqrt(2);0;1/sqrt(2)];
cos_al = 0.9;
alpha = acos(cos_al);
Rpts = randpts(Nr,alpha,0);
c = [ct(1);ct(2);ct(3)];
vc = [-c(2);c(1);0]/sqrt(c(1)^2+c(2)^2);
vcc = cross(vc,c);
M = [vcc vc c]; % rotation matrix
XR = M*Rpts';      % now XR is centered at c

% compute Fourier coefficients
load DH_64; % Driscoll-Healy quadrature X_W of 2^{14} points
X = X_W(:,1:3);
f = (X(:,1).*cos(X(:,3)) - X(:,3).*sin(X(:,3)) + 0.5*cos(X(:,3))).*exp(X(:,1));
FH = Fhat(X_W,f',50);
goodU32 = cap_spectral(FH,32,1,XR);
exact_sol = exp(XR(1,:)).* cos(XR(3,:));
err = abs(goodU32-exact_sol); % which is 1.199e-14 by Matlab 7.0 at UNSW.
