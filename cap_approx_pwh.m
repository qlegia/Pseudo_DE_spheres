function [err,err_noh,exact_err,g_exact_err] = cap_approx_pwh(goodFH,goodL,FH,Lsol)
% Input:  
%         FH =  widehat{f}_{ell,k} produced by new_spectral
%         Lsol = max degree of the spherical harmonics
%         goodFH = widehat{f}_{ell,k} of the good solution
%         goodL = max degree of sph. harm. in good solution
%              we assume goodL >= Lsol
%         SolPts : (3, N) matrix of Cartesian coordinates 
%           on which we want to compute the solution    
% Output:
%         local solution on the set of points SolPts
%         soln is a (1,N) matrix
% Purpose: 
%         measure local error on spherical caps
% -------------------------------------------------------
% Q.T. Le Gia, UNSW, Sydney, Australia.
% 30-Aug-2005.
% -------------------------------------------------------
% The eigenvalues of the pseudo-differential equation
lambda = ones(1,goodL+1); % identity operator
max_q = 7;
err = zeros(1,max_q);
g_hvals = zeros(max_q,goodL+1);
g_samp = [0:goodL]./goodL;
samp = [0:Lsol]./Lsol;
hvals = zeros(max_q,Lsol+1);
for q=1:max_q
  g_hvals(q,:) = powhfunc(q,g_samp);
  hvals(q,:) = powhfunc(q,samp);
  g_hvals(q,1) = 1; 
  hvals(q,1) = 1;
end
% set the eigenvalues of the pseudo-diff operator here
%lambda(1) = 0.5; %differential operator
%lambad(1) = 2.0; %integral operator
%for ell=1:Lsol
    %lambda(ell+1) = (ell+0.5); %differential operator
    %lambda(ell+1) = (2*ell+2)/(2*ell+1); % integral operator
%end     
% load pre-compute sph. harmonics
load Y0_100_XR2;
% begin to compute solutions
goodFHY  = cell(q,1);  
FHY = cell(q,1);
for ell=0:goodL
  r_indx = 1 + ell*(ell+1);
  Y = Y0_100(r_indx-ell:r_indx+ell,:);
  if (ell <= Lsol)
   for m=-ell:ell
      g_FHYnoh = goodFH(r_indx+m)*Y(m+ell+1,:)/lambda(ell+1);
      FHYnoh = FH(r_indx+m)*Y(m+ell+1,:)/lambda(ell+1);
      for q=1:max_q
        goodFHY{q}(r_indx+m,:) = g_hvals(q,ell+1)*goodFH(r_indx+m)*Y(m+ell+1,:)/lambda(ell+1);
        FHY{q}(r_indx+m,:) = hvals(q,ell+1)*FH(r_indx+m)*Y(m+ell+1,:)/lambda(ell+1);
      end  
   end
  else % only worry about good solution
   for m=-ell:ell
      g_FHYnoh = goodFH(r_indx+m)*Y(m+ell+1,:)/lambda(ell+1);
      for q=1:max_q
        goodFHY{q}(r_indx+m,:) = g_hvals(q,ell+1)*goodFH(r_indx+m)*Y(m+ell+1,:)/lambda(ell+1);
      end
   end
  end %endif
end
exactf = non_smooth_f(XR2');
exactf = exactf';
% solutions without function h
soln_noh = sum(FHYnoh);
goodsoln_noh = sum(g_FHYnoh);
err_noh = max(abs(soln_noh-goodsoln_noh));
% compute errors
for q=1:max_q
   soln(q,:) = sum(FHY{q}); % sum up all ell's and m's at each point for the solution
   goodsol(q,:) = sum(goodFHY{q});
   err(q) = max(abs(soln(q,:)-goodsol(q,:)));
   exact_err(q) = max(abs(goodsol(q,:)-exactf));
   g_exact_err(q) = max(abs(soln(q,:)-exactf));
end  

