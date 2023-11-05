function soln = cap_spectral(FH,Lsol,with_h,SolPts)
% Input:  
%         FH =  widehat{f}_{ell,k} produced by new_spectral
%         Lsol = max degree of the spherical harmonics
%         with_h = 0  no mask function h
%         with_h > 0  mask fn h with smoothness with_h   
%         SolPts : (3, N) matrix of Cartesian coordinates 
%           on which we want to compute the solution    
% Output:
%         local solution on the set of points SolPts
%         soln is a (1,N) matrix
% Purpose: 
%         measure local error on spherical caps
% -------------------------------------------------------
% Q.T. Le Gia, UNSW, Sydney, Australia.
% 29-Aug-2005.
% -------------------------------------------------------
% The eigenvalues of the pseudo-differential equation
lambda = zeros(1,Lsol+1);
hvals = zeros(1,Lsol+1);
samp = [0:Lsol]./Lsol;
if (with_h>0)
  hvals = hfunc(with_h,samp);
  hvals(1) = 1; 
end
lambda(1) = 0.5;
for ell=1:Lsol
    lambda(ell+1) = (ell+0.5);
end     
% begin to compute solutions
for ell=0:Lsol
  [theta,phi] = x2theta(SolPts); 
  Y = hom_harmonics(ell,theta,phi,'real');
  r_indx = 1 + ell*(ell+1);
  if (with_h > 0)
    for m=-ell:ell
      FHY(r_indx+m,:) = hvals(ell+1)*FH(r_indx+m)*Y(m+ell+1,:)/lambda(ell+1);
    end
  else
    for m=-ell:ell
      FHY(r_indx+m,:) = FH(r_indx+m)*Y(m+ell+1,:)/lambda(ell+1);
    end
  end
end
soln = sum(FHY); % sum up all ell's and m's at each point for the solution

