function [err,err_noh] = cap_approx_exact(Y,FH,Lsol,SolPts,exactf)
% Input:  
%         FH =  widehat{f}_{ell,k} produced by new_spectral
%         Lsol = max degree of the spherical harmonics         
%         SolPts : (3, N) matrix of Cartesian coordinates 
%           on which we want to compute the solution    
% Output:
%         err of local solution on the set of points SolPts
%         soln is a (1,N) matrix
% Purpose: 
%         measure local error on spherical caps
% -------------------------------------------------------
% Q.T. Le Gia, UNSW, Sydney, Australia.
% 30-Aug-2005.
% -------------------------------------------------------
% The eigenvalues of the pseudo-differential equation
lambda = ones(1,Lsol+1); % identity operator
max_q = 7;
err = zeros(1,max_q);
samp = [0:Lsol]./Lsol;
hvals = zeros(max_q,Lsol+1);
for q=1:max_q
  hvals(q,:) = hfunc(q,samp); 
  hvals(q,1) = 1;
end
% set the eigenvalues of the pseudo-diff operator here
%lambda(1) = 0.5; %differential operator
%lambad(1) = 2.0; %integral operator
%for ell=1:Lsol
    %lambda(ell+1) = (ell+0.5); %differential operator
    %lambda(ell+1) = (2*ell+2)/(2*ell+1); % integral operator
%end     
% begin to compute solutions
FHY = cell(q,1);
for ell=0:Lsol
  r_indx = 1 + ell*(ell+1);
  for m=-ell:ell
      for q=1:max_q
        FHY{q}(r_indx+m,:) = hvals(q,ell+1)*FH(r_indx+m)*Y(r_indx+m,:)/lambda(ell+1);
      end  
      FHYnoh(r_indx+m,:) = FH(r_indx+m)*Y(r_indx+m,:)/lambda(ell+1);
   end
end
soln = zeros(q,length(SolPts));
for q=1:max_q
   soln(q,:) = sum(FHY{q}); % sum up all ell's and m's at each point for the solution
   err(q) = max(abs(soln(q,:)-(exactf(:))'));
end  
soln_noh = sum(FHYnoh);
err_noh = max(abs(soln_noh-(exactf(:))'));
