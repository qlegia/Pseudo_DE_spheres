function FH = Fhat(X_W,f,L)
% Input  : the quadrature X_W, the function f, 
%          L : all ell=0..L-1 will be considered
% Output : the Fourier coefficients widehat{f}_{ell,m}
%          FH(1+ell*(ell+1)+m) = widehat{f}_{ell,m}
% Purpose: compute widehat{f}_{ell,k} using the quadrature X_W
% ------------------------------------------------------------------------
% Q.T. Le Gia, UNSW, 29-Aug-2005.
% ------------------------------------------------------------------------
% convert X_W to spherical coordinates of points and extract the weights
X = (X_W(:,1:3))'; % now X is a 3 by N matrix
N = length(X);
W = X_W(:,4)'; % transpose so that the weights W is a row vector
[theta,phi] = x2theta(X);
tic
FH = zeros(L^2,1); % zeros column vector
for ell=0:L-1
   % homogeneous spherical harmonics of order ell
   ell
   Y = hom_harmonics(ell,theta,phi,'real');
   r_indx = 1 + ell*(ell+1);
   for m=-ell:ell
      % compute the Fourier coefficients fhat(ell,m) by quadrature then 
      % store into FH(r_indx+m) 
      FH(r_indx+m) = sum(Y(m+ell+1,:) .* f .* W);        
   end
   toc
end
