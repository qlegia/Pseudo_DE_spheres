% Test the Neumann problem with Driscoll-Healy quadrature rule
% This code use to generate the error for Table 2 in mhaskar_report.tex
% Quoc Thong Le Gia, UNSW, Sydney, Australia.
% Oct-28-2004
% -----------------------------------------------------------------------------
clear all; 
%load DH_16;  %  1024 points
%load DH_32;  %  4096 points
%load DH_64;  % 16384 points
load DH_128; % 65536 points
X = X_W(:,1:3);
f = 2*ones(length(X),1);
[FH,soln] = integral2(X_W,f',42,40,0);
[Sx,Sy,Sz] = sphere(50);
err = soln - 1;
max(max(abs(err)))
norm(abs(err))/51
[FH,soln] = integral2(X_W,f',42,40,1);
[Sx,Sy,Sz] = sphere(50);
err = soln - 1;
max(max(abs(err)))
norm(abs(err))/51