function [err,err_noh] = cap_int_pwh1(goodsoln,goodsoln_noh,FH,Lsol)
% Input:  
%         FH =  widehat{f}_{ell,k} produced by Fhat.m 
%         Lsol = max degree of the spherical harmonics 
%         goodsoln = pre-computed good solution of high approx. degree
%         SolPts : (3, N) matrix of Cartesian coordinates 
%           on which we want to compute the solution    
% Output:
%         local solution on the set of points SolPts
%         soln is a (1,N) matrix
% Purpose: 
%         use powhfunc.m
%         measure local error on spherical caps
% -------------------------------------------------------
% Q.T. Le Gia, UNSW, Sydney, Australia.
% 2-Sep-2005.
% -------------------------------------------------------
% The eigenvalues of the pseudo-differential equation
lambda = ones(1,Lsol+1); % identity operator
max_q = 7;
err = zeros(1,max_q);
samp = [0:Lsol]./Lsol;
hvals = zeros(max_q,Lsol+1);
for q=1:max_q
  hvals(q,:) = powhfunc(q,samp);
  hvals(q,1) = 1;
end
% set the eigenvalues of the pseudo-diff operator here
%lambda(1) = 0.5; %differential operator
lambda(1) = 2.0; %integral operator
for ell=1:Lsol
    %lambda(ell+1) = (ell+0.5); %differential operator
    lambda(ell+1) = (2*ell+2)/(2*ell+1); % integral operator
end     
% load pre-compute sph. harmonics from 0 to 100
load Y0_100_XR2;
% begin to compute solutions
FHY = cell(q,1);
for ell=0:Lsol
  r_indx = 1 + ell*(ell+1);
  Y = Y0_100(r_indx-ell:r_indx+ell,:);
  for m=-ell:ell
      FHYnoh = FH(r_indx+m)*Y(m+ell+1,:)/lambda(ell+1);
      for q=1:max_q
        FHY{q}(r_indx+m,:) = hvals(q,ell+1)*FH(r_indx+m)*Y(m+ell+1,:)/lambda(ell+1);
      end  
  end
end
% solutions without function h
soln_noh = sum(FHYnoh);
err_noh = max(abs(soln_noh-goodsoln_noh));
% compute errors
for q=1:max_q
   soln(q,:) = sum(FHY{q}); % sum up all ell's and m's at each point for the solution
   err(q) = max(abs(soln(q,:)-goodsoln(q,:)));
end  
