%------------------------------------------------------------------
% Generate require errors for pseudo-spectral method.
% Assume that we know all the eigenvalues lamda_ell of
%    the pseudo-differential operators
% The Fourier coefficients has been computed using FFT_Fhat.m
% -----------------------------------------------------------------
% Q.T. Le Gia, UNSW, Sydney, Australia.
% Sep-13-2005
% -----------------------------------------------------------------
function [sol,sol_noh] = FFT_pseudo_fn(Lsol,cfh,XR,equ_type)
% cfh contains complex Fourier coefficients of hat{f}(ell,m)
% for all ell <= Lsol
% Lsol: approximate solution will be a combination of all sph
% harm of deg < Lsol 
l_range = (0:Lsol);

% set the eigenvalues of the pseudo-operators here ...
lambda = ones(1,Lsol+1); % identity operator by default
if (equ_type == 1)
  % differential operator
  lambda = (0:Lsol)+1/2;
else if (equ_type == 2)
  % Green integral operator
  lambda = 1./((0:Lsol)+0.5).^2;
else if (equ_type == 3)
  % Integral equ of the 2nd kind
  lambda = (2*l_range+1)./(2*l_range+2);
end

tic
max_q=5;
% all the spherical harmonics at given points
if (size(XR,2)<=3)
    XR = XR';
end    
Nr = length(XR);
[theta,phi] = x2theta(XR);
sol_m = zeros(Lsol,Nr);
sol = zeros(max_q,Nr);
sol_m_noh = zeros(Lsol,Nr);
sol_noh = zeros(max_q,Nr);
sol00 = cfh(1)/sqrt(4*pi)/lambda(1);
for q=1:max_q
   hsol = hfunc2(q,(0:Lsol-1)./Lsol);
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
  sol(q,:) = sol(q,:) + hsol(1)*sol00; 
  % add the constant sol00 = u_{0,0} to the solution
  sol_noh = sol_noh + sol00;
  ell
  toc
end
