% test_cap_spectral
% 29-Aug-2005
% Q.T. Le Gia, UNSW, Sydney, Australia.
clear all;
tic
% we want solution on a spherical cap centered at ctradius alpha
Nr = 10000; %number of random points for ssolutions 
with_h = 5; % smoothness of the function h

c1 = [1/sqrt(2);0;1/sqrt(2)];            cos_al1 = 0.9;
c2 = [-1/sqrt(2);0;1/-sqrt(2)];          cos_al2 = 0.9;
c3 = c2;                                 cos_al3 = 0.93;
c4 = [-1/sqrt(3);-1/sqrt(3);-1/sqrt(3)]; cos_al4 = 0.9;

alpha1 = acos(cos_al1);
alpha2 = acos(cos_al2);
alpha3 = acos(cos_al3);
alpha4 = acos(cos_al4);

XR1 = ranpts_cap(Nr,c1,alpha1);
XR2 = ranpts_cap(Nr,c2,alpha2);
XR3 = ranpts_cap(Nr,c3,alpha3);
XR4 = ranpts_cap(Nr,c4,alpha4);
% compute Fourier coefficients for the good solution
load DH_64; % Driscoll-Healy quadrature X_W of 2^{14} points
X = X_W(:,1:3);
% f = (X(:,1).*cos(X(:,3)) - X(:,3).*sin(X(:,3)) + 0.5*cos(X(:,3))).*exp(X(:,1));
f = non_smooth_f(X);
goodFH = Fhat(X_W,f',50);
% compute Fourier coeffs for the approximation
clear X_W;
clear X;
clear f;
load DH_32;
X = X_W(:,1:3);
%f = (X(:,1).*cos(X(:,3)) - X(:,3).*sin(X(:,3)) + 0.5*cos(X(:,3))).*exp(X(:,1));
f = non_smooth_f(X);
FH = Fhat(X_W,f',27);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for with_h = 0:1
  % compute the good solution
  goodU32 = zeros(4,Nr);
  goodU32(1,:) = cap_spectral(goodFH,32,with_h,XR1);
  goodU32(2,:) = cap_spectral(goodFH,32,with_h,XR2);
  goodU32(3,:) = cap_spectral(goodFH,32,with_h,XR3);
  goodU32(4,:) = cap_spectral(goodFH,32,with_h,XR4);
  % compute the approximation
  soln = zeros(4,Nr);
  soln(1,:) = cap_spectral(FH,13,with_h,XR1);
  soln(2,:) = cap_spectral(FH,13,with_h,XR2);
  soln(3,:) = cap_spectral(FH,13,with_h,XR3);
  soln(4,:) = cap_spectral(FH,13,with_h,XR4);
  % now the errors for each spherical cap
  err = zeros(1,4);
  err(1) = max(abs(goodU32(1,:)-soln(1,:)));
  err(2) = max(abs(goodU32(2,:)-soln(2,:)));
  err(3) = max(abs(goodU32(3,:)-soln(3,:)));
  err(4) = max(abs(goodU32(4,:)-soln(4,:)));
  sprintf('%0.10e %0.10e %0.10e %0.10e',err(1),err(2),err(3),err(4))
end  
toc