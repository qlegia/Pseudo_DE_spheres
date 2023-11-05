function XR = ranpts_cap(N,ct,r)
%Written by Dr. Quoc Thong Le Gia on August 29, 2005
%call XR=randpts_cap(N,ct,alpha)
% generate N random points on a spherical caps centered ct, radius r
%The cap is given by $\{x: x\cdot ct\ge r\}$
%r must be in the range [-1,1]
%uses randpts.m

%First generate N random points on the cap of radius r, centered at the north pole
alpha=acos(r);
Rpts=randpts(N,alpha,0);
%Now rotate these points to get centered at ct
c = [ct(1);ct(2);ct(3)];
vc = [-c(2);c(1);0]/sqrt(c(1)^2+c(2)^2);
vcc = cross(vc,c);
M = [vcc vc c]; % rotation matrix
XR = M*Rpts';   
