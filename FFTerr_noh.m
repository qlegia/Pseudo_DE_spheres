% Sep-9-2005
% 
% Generate require errors for pseudo-spectral method.
% Assume that we know all the eigenvalues lamda_ell of
%    the pseudo-differential operators
% The Fourier coefficients has been computed using FFT_Fhat.m
% --------------------------------------------------------------
% Q.T. Le Gia, UNSW, Sydney, Australia.
% --------------------------------------------------------------
Lmax = 1024; % the ideal solution
% set the eigenvalues of the pseudo-operators here ...
% differential operator
%lambda = (0:Lmax)+1/2;
%integral operator
lambda = 1./((0:Lmax)+1/2).^2;
% load pre-computed Fourier coefficients
load good_cfh; 
Lsol = 256; % the max degree of the approximate solution
cfh = cfh256;
% the mask
max_q=5;
idx = 1:Lsol^2;
idx_t = Lsol^2+1:Lmax^2;
for q=1:max_q
   hsol = powhfunc(q,(0:Lsol-1)./Lsol);
   hLmax = powhfunc(q,(0:Lmax-1)./Lmax);
   for ell=0:Lmax-1
      lower_b = ell*(ell+1)+1-ell;
      upper_b = ell*(ell+1)+1+ell;
      harr_Lmax(lower_b:upper_b) = hLmax(ell+1);
      arr_lambda(lower_b:upper_b) = lambda(ell+1);
      if (ell< Lsol)
           hsol = powhfunc(q,(0:Lsol-1)./Lsol);
           harr_Lsol(lower_b:upper_b) = hsol(ell+1);
      end % endif
   end
  err_head(idx) = (harr_Lmax(idx).*cfh1024(idx) - harr_Lsol(idx).*cfh(idx))./arr_lambda(idx);
  err_tail(1:length(idx_t)) =  harr_Lmax(idx_t).*cfh1024(idx_t)./arr_lambda(idx_t); 
  err1(q) = sum(err_head)
  err2(q) = sum(err_tail)
  err(q) = err1(q) + err2(q)
end
