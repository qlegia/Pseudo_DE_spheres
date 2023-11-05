function [goodsol,goodsoln_noh] = cap_diff_pwh(goodFH,goodL)
% Input:  
%         goodFH = widehat{f}_{ell,k} of the good solution
%         goodL = max degree of sph. harm. in good solution
%         SolPts : (3, N) matrix of Cartesian coordinates 
%           on which we want to compute the solution    
% Output:
%         local solution on the set of points SolPts
%         soln is a (1,N) matrix
% Purpose: 
%         measure local error on spherical caps
% -------------------------------------------------------
% Q.T. Le Gia, UNSW, Sydney, Australia.
% 2-Sep-2005.
% -------------------------------------------------------
% The eigenvalues of the pseudo-differential equation
lambda = ones(1,goodL+1); % identity operator
max_q = 5;
err = zeros(1,max_q);
g_samp = (0:(goodL-1))./goodL;
for q=1:max_q
  g_hvals(q,:) = powhfunc(q,g_samp);
  g_hvals(q,1) = 1; 
end
% set the eigenvalues of the pseudo-diff operator here
lambda(1) = 0.5; %differential operator
%lambad(1) = 2.0; %integral operator
for ell=1:goodL-1
    lambda(ell+1) = (ell+0.5); %differential operator
    %lambda(ell+1) = (2*ell+2)/(2*ell+1); % integral operator
end     
% load pre-compute sph. harmonics
load Y0_100_XR2;
tic
% begin to compute solutions
goodFHY  = cell(q,1);  
for ell=0:goodL-1
  if (ell>60)
     ell
     toc
     if (ell>100)
         clear Y0_100_XR2;
         load Y101_150_XR2;
      end
  end
  r_indx = 1 + ell*(ell+1);
  if (ell<=100)
     Y = Y0_100(r_indx-ell:r_indx+ell,:);
  else
     r_indx = r_indx-10201; % (100+1)^2=10201
     toc
     Y = Y101_150(r_indx-ell:r_indx+ell,:);
  end

  for m=-ell:ell
      g_FHYnoh(r_indx+m,:) = goodFH(r_indx+m)*Y(m+ell+1,:)/lambda(ell+1);
      for q=1:max_q
        goodFHY{q}(r_indx+m,:) = g_hvals(q,ell+1)*goodFH(r_indx+m)*Y(m+ell+1,:)/lambda(ell+1);
      end  
  end
end
% solutions without function h
goodsoln_noh = sum(g_FHYnoh);
% compute errors
for q=1:max_q
   goodsol(q,:) = sum(goodFHY{q});
end 
