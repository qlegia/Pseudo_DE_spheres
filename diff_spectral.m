function [soln] = diff_spectral(FH,Lsol,with_h)
% diff_spectral(X_W,f,Lsol,with_h)
% --------------------------------------------------------------------------
% spectral method using Fourier transform
% f is a row vector, X_W contains the quadrature rule as Nx4 matrix
% Lsol is the highest degree of spherical harmonics in the solution
% with_h > 0 : turn on the filter function h
% with_h = 0 : turn off the filter function h
% --------------------------------------------------------------------------
% Quoc Thong Le Gia, University of N.S.W., Sydney, Australia. 
% 04-Oct-2004 : first draft
% 28-Aug-2005 : bug on mask function hval(x) fixed
% 17-Sep-2005 : more vectorized 
% --------------------------------------------------------------------------
%
FH = FH(1:Lsol^2); % row vector
% The eigenvalues of the pseudo-differential equation
lambda = zeros(1,Lsol^2);
if (with_h>0)
  hLsol = hfunc2(with_h,[0:1/Lsol:1]);
end
hvals = zeros(1,Lsol^2); 
lambda(1) = 0.5;
for ell=1:Lsol-1
     r_index = 1+ell*(ell+1);
     lambda(r_index-ell:r_index+ell) = ones(1,2*ell+1)*(ell+0.5);
end
if (with_h>0)
  hvals(1) = 1;
  for ell=1:Lsol-1  
     r_index=1+ell*(ell+1); 
     hvals(r_index-ell:r_index+ell) = ones(1,2*ell+1)*hLsol(ell);
  end   
end     
% we want a good looking solution
num = 50;
[Sx,Sy,Sz] = sphere(num);
soln = zeros(num+1,num+1);
%  Ysol will be a grid corresponding to Sx,Sy,Sz 
for k=1:num+1
    [theta,phi]=x2theta([Sx(k,:);Sy(k,:);Sz(k,:)]);
    Ysol = s_harmonics(Lsol-1,theta,phi,'real');
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
