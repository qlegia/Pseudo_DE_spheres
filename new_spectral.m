function [FH,soln] = new_spectral(X_W,f,L,Lsol,with_h)
% new_spectral(X_W,f,L,Lsol,with_h)
% --------------------------------------------------------------------------
% spectral method using Fourier transform
% f is a row vector, X_W contains the quadrature rule as Nx4 matrix
% L is highest order of spherical harmonics to be computed
% Lsol is the highest degree of spherical harmonics in the solution
% with_h > 0 : turn on the filter function h
% with_h = 0 : turn off the filter function h
% --------------------------------------------------------------------------
% Quoc Thong Le Gia, University of N.S.W., Sydney, Australia. 
% 04-Oct-2004 : first draft
% 28-Aug-2005   : bug on mask function hval(x) fixed 
% 29-Aug-2005 : hvals(1) = 1
% --------------------------------------------------------------------------
% convert X_W to spherical coordinates of points and extract the weights
X = (X_W(:,1:3))'; % now X is a 3 by N matrix
N = length(X);
W = X_W(:,4)'; % transpose so that the weights W is a row vector
[theta,phi] = x2theta(X);
FH = zeros(L^2,1); % zeros column vector
for ell=0:L-1
   % homogeneous spherical harmonics of order ell
   Y = hom_harmonics(ell,theta,phi,'real');
   r_indx = 1 + ell*(ell+1);
   for m=-ell:ell
      % compute the Fourier coefficients fhat(ell,m) then store into FH(r_indx+m)
      FH(r_indx+m) = sum(Y(m+ell+1,:) .* f .* W);  % using the quadrature rule      
   end  
end  
% The eigenvalues of the pseudo-differential equation
lambda = zeros(1,(Lsol+1)^2);
samp = [1:Lsol]./Lsol;
hLsol = hfunc(with_h,samp);
hvals = zeros(1,(Lsol+1)^2);
hvals(1) = 1; 
lambda(1) = 0.5;
for ell=1:Lsol
     r_index = 1+ell*(ell+1);
     lambda(r_index-ell:r_index+ell) = ones(1,2*ell+1)*(ell+0.5);
     hvals(r_index-ell:r_index+ell) = ones(1,2*ell+1)*hLsol(ell);
end     
% we want a good looking solution
num = 50;
[Sx,Sy,Sz] = sphere(num);
soln = zeros(num+1,num+1);
%  Ysol will be a grid corresponding to Sx,Sy,Sz 
for k=1:num+1
    [theta,phi]=x2theta([Sx(k,:);Sy(k,:);Sz(k,:)]);
    Ysol = s_harmonics(Lsol,theta,phi,'real'); % compute all the spherical harmonics
    d_sol = (Lsol+1)^2;                        % the number of rows of matrix Ysol
    FH3 = zeros((Lsol+1)^2,num+1);
    if (with_h > 0) 
        for n=1:d_sol
            FH3(n,:) = hvals(n)*FH(n)*Ysol(n,:)/lambda(n);      
        end
    else
        for n=1:d_sol
            FH3(n,:) = FH(n)*Ysol(n,:)/lambda(n);      
        end
    end    
    for j=1:num+1
        soln(k,j) = sum(FH3(:,j));    
    end    
end    
