function [FH,soln] = MCspectral(X,f,Lsol,with_h)
% MCspectral(X_W,f,Lsol,with_h)
% --------------------------------------------------------------------------
% spectral method using Fourier transform and Monte-Carlo method
% f is a row vector, X_W contains the quadrature rule as Nx4 matrix
% Lsol is the highest degree of spherical harmonics in the solution
% with_h = 1 : turn on the filter function h
% with_h = 0 : turn off the filter function h
% --------------------------------------------------------------------------
% Quoc Thong Le Gia, University of NSW, Sydney, Australia, 22-Aug-2005.
% Sep-16-2005: vectorize FH(r_indx+m) to FH(r_indx+(-ell:ell))
%              vectorize FH3
% --------------------------------------------------------------------------
% convert X to spherical coordinates of points and extract the weights
X = X'; % now X is a 3 by N matrix
N = length(X);
if (N<1000)
  sprintf('ERROR')
  return
end
[theta,phi] = x2theta(X);
FH = zeros(1,Lsol^2); % zeros column vector
for ell=0:Lsol-1
   % homogeneous spherical harmonics of order ell
   Y = hom_harmonics(ell,theta,phi,'real');
   r_indx = 1 + ell*(ell+1);
   m_range = (-ell:ell);
   FH(r_indx+m_range) = monte_carlo1(N,Y,f);
end  
% The eigenvalues of the pseudo-differential equation
lambda = zeros(1,Lsol^2);
lambda(1) = 0.5;
hval = hfunc2(1,[0:1/Lsol:(Lsol-1)]);
hval(1) = 1; % just make sure
for ell=1:Lsol-1
     r_index = 1+ell*(ell+1);
     lambda(r_index-ell:r_index+ell) = ones(1,2*ell+1)*(ell+0.5);
     hvals(r_index-ell:r_index+ell) = ones(1,2*ell+1)*hval(ell+1);     
end     
% we want a good looking solution
num = 50;
[Sx,Sy,Sz] = sphere(num);
soln = zeros(num+1,num+1);
%  Ysol will be a grid corresponding to Sx,Sy,Sz 
for k=1:num+1
    [theta,phi]=x2theta([Sx(k,:);Sy(k,:);Sz(k,:)]);
    Ysol = s_harmonics(Lsol-1,theta,phi,'real');
     % compute all the spherical harmonics
    FH3 = zeros(Lsol^2,num+1);
    if (with_h == 1) 
        FHld = diag(hvals.*FH./lambda);
        FH3 = FHld*Ysol;
    else
        FHld = diag(FH./lambda);
        FH3 = FHld*Ysol;      
    end    
    soln(k,:) = sum(FH3);
end    
