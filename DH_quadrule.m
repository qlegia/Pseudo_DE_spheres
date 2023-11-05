function X_W = DH_quadrule(k)
% DH_quadrule(k) computes Driscoll-Healy quadrature rule on the sphere S^2 
% for n=2^k samples along the z axis. The quadrature rule gives exact Fourier 
% coefficients for spherical harmonics up to degree n/2. 
% The result X_W is a N=n^2 by 4 matrix contains Cartesian points and weights.
% ---------------------------------------------------------------------------
% Quoc Thong Le Gia, Texas A & M University, 6/15/2004.
% ---------------------------------------------------------------------------
n = 2^k;
W = zeros(1,n);
for j=0:n-1    
  v_ell = 2*[0:n/2-1]+1;    
  W(j+1) = (2*sqrt(2)/n)*sin(pi*j/n)*sum(sin(v_ell*pi*j/n)./v_ell);
end
  theta = (pi/n)*[0:1:n-1]';
  phi = (2*pi/n)*[0:1:n-1];
  sintheta = sin(theta);  
  X = sintheta*cos(phi);
  Y = sintheta*sin(phi);
  Z = cos(theta)*ones(1,n);
  W = W*(sqrt(2*pi)/n)*2*sqrt(pi); 
  % adjust normalized constant 2*sqrt(pi)
  X_W = [X(:,1) Y(:,1) Z(:,1) W'];
for j=2:n    
  X_W = [X_W; X(:,j) Y(:,j) Z(:,j) W'];
end    
