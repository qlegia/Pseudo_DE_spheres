% Sep-9-2005
% Sep-10-2005: added line 57, adjust the constant 
%   lambda^{-1}_0 hat{f}_{0,0} Y_{0,0}
%------------------------------------------------------------------
% Generate require errors for pseudo-spectral method.
% Assume that we know all the eigenvalues lamda_ell of
%    the pseudo-differential operators
% The Fourier coefficients has been computed using FFT_Fhat.m
% -----------------------------------------------------------------
% Q.T. Le Gia, UNSW, Sydney, Australia.
% -----------------------------------------------------------------
clear all;
format long;
Lsol = 256;
% set the eigenvalues of the pseudo-operators here ...
% differential operator
%lambda = (0:Lmax)+1/2;
%integral operator
lambda = 1./((0:Lsol)+0.5).^2;
% lambda = ones(1,Lsol+1); %id operator
% load pre-computed Fourier coefficients
%load cfh1024_256;
load cfh2048_256;
% pre-computed hat{f}_{ell,m} using quadrature exact upto degree 255
% 
tic
max_q=5;
% all the spherical harmonics at given points
Nr = 1000; % 1000 random points on a cap
%XR = ranpts_cap(1000,[-1/sqrt(2);0;-1/sqrt(2)],acos(0.9));
load XR2;
if (size(XR2,2)<=3)
    XR2 = XR2';
end    
XR = XR2;
[theta,phi] = x2theta(XR);
sol_m = zeros(Lsol,Nr);
sol = zeros(max_q,Nr);
sol_m_noh = zeros(Lsol,Nr);
sol_noh = zeros(max_q,Nr);
sol00 = cfh(1)/sqrt(4*pi)/lambda(1);
for q=1:max_q
   hsol = powhfunc(q,(0:Lsol-1)./Lsol);
   % treat ell=0, m=0 separately
   % sol(q,1:Nr) = (cfh(1)/sqrt(4*pi)/lambda(1))*ones(1,Nr);
   for ell=1:Lsol-1
      lower_r = ell*(ell+1)+1-ell;
      center_r = ell*(ell+1)+1;
      upper_r = ell*(ell+1)+1+ell;
      Pell = legendre(ell,XR(3,:),'sch');
      pl_minus = (-1).^(1:ell);
      phase = exp(i*(1:ell)'*phi)./sqrt(2);
      yy_plus  = Pell(2:ell+1,:).*phase;
      yy_0 = Pell(1,:);
      yy_minus = flipud(conj(yy_plus));
      Y = (sqrt(2*ell+1)/sqrt(4*pi)).* [yy_minus;yy_0;diag(pl_minus)*yy_plus];     
      sYell = cfh(lower_r:upper_r)*Y(1:2*ell+1,:);
      sol_m(ell+1,:) = (hsol(ell+1)/lambda(ell+1)).*sYell;
      sol_m_noh(ell+1,:) = sYell./lambda(ell+1);
      %sol_m(ell+1,:) = sYell;   
  end
  sol(q,1:Nr) = sum(sol_m); % sum over all ell
  sol_noh = sum(sol_m_noh);   
  % add the constant sol00 = u_{0,0}
  sol(1,:) = sol(1,:) + hsol(1)*sol00; 
  % add the constant sol00 = u_{0,0} to the solution
  sol_noh = sol_noh + sol00;
  ell
  toc
end
% now compute the errors
load ideal_sol;
for q=1:max_q
    err(q) = max(abs(ideal-real(sol(q,:))));
end  
err_noh = max(abs(ideal-real(sol_noh)));
display(err')
display((err-sol00)')
display(err_noh)