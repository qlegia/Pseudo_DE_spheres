function max_e = cap_err(N,err,w,cos_alpha)
% Input  : err is N by N matrix representing a grid of values on the 
%          surface of the sphere,
%          w is a (3,1) vector gives the center of the cap, 
%          alpha is the spherical radius of the spherical cap   
% Output : actual maximum error on the cap
% Purpose: calculate the error on a local spherical caps
% ----------------------------------------------------------------------------
% Q.T.Le Gia, UNSW, Aug-10-2005.
% ----------------------------------------------------------------------------
% generating the spherical grid
[Sx,Sy,Sz] = sphere(N-1);
ctheta = zeros(N,N);
cap_e = zeros(1,N);
for k = 1:N
   ctheta(1:N) = w(1)*Sx(k,:)+w(2)*Sy(k,:)+w(3)*Sz(k,:);
   idx = find(ctheta>cos_alpha);
   if (length(idx)>0) 
     cap_e(k) = max(err(k,idx));
   end
end
max_e = max(cap_e);
