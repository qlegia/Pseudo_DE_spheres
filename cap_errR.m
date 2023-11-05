function max_cap = cap_errR(err,w,cos_alpha)
% Input:    err is (N,3) matrix representing errors on the sphere
%           w is a (3,1) vector indicates center of the spherical cap
%           alpha is the center of the spherical cap
% Output:   maximum error on the cap
% Purpose:  calculate error on local spherical caps using random points
%           on the cap centered at w radius alpha
%------------------------------------------------------------------------
% Q.T.Le Gia, UNSW, 29-Aug-2005.
%------------------------------------------------------------------------
N = 10000;
%generate N random points on a cap center at north pole
XR = randpts(N);
% rotate those N points to the center w
X = M*XR';

