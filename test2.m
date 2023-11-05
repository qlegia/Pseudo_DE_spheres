% Test the Neumann problem with Driscoll-Healy quadrature rule
% This code use to generate the error for Table 2 in mhaskar_report.tex
% Quoc Thong Le Gia, UNSW, Sydney, Australia.
% Oct-5-2004
% -----------------------------------------------------------------------------
clear all; 
%load DH_16;  %  1024 points
%load DH_32;  %  4096 points
%load DH_64;  % 16384 points
load DH_128; % 65536 points
X = X_W(:,1:3);
f = ( X(:,1).*cos(X(:,3)) - X(:,3).*sin(X(:,3))+ 0.5*cos(X(:,3)) ).*exp(X(:,1));
[FH,soln] = new_spectral(X_W,f',42,40,0);
% u_exact = exp(X(:,1)).*cos(X(:,3));
[Sx,Sy,Sz] = sphere(50);
err = soln - exp(Sx).*cos(Sz);
max(max(abs(err)))
norm(abs(err))/51
[FH,soln] = new_spectral(X_W,f',42,40,1);
% u_exact = exp(X(:,1)).*cos(X(:,3));
[Sx,Sy,Sz] = sphere(50);
err = soln - exp(Sx).*cos(Sz);
max(max(abs(err)))
norm(abs(err))/51
