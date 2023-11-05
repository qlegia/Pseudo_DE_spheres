% MATLAB code used to generate Table 1 in mhaskar_report.tex
% Q. T. Le Gia, Oct-6-2004.
% Sep-17-2005: modified with diff_spectral.m
% ------------------------------------------------------------------------------
clear all;
%-------------------------------------------------------------------------------
%f = (X(:,1).*cos(X(:,3)) - X(:,3).*sin(X(:,3)) + 0.5*cos(X(:,3))).*exp(X(:,1));
FH = realFhat_fs(8,64);
sol = diff_spectral(FH,64,2);
% u_exact = exp(X(:,1)).*cos(X(:,3));
% now we have to compute the error between the exact solution and 
% the approximation
num = 50;
[Sx,Sy,Sz] = sphere(num);
err_max = zeros(num+1,num+1);
for k=1:num+1       
       err_max(k,:) = abs(sol(k,:) - exp(Sx(k,:)).*cos(Sz(k,:)));               
end    
max(max(err_max))
norm(err_max)/(num+1)
%-------------------------------------------------------------------------------
% Results of Sep-17-2005
% Number of D-H quadrature points:  (128)^2=16384 points
%   err_noh = max|U32 - u_exact| = 1.498801083243961e-14
% ----------------------------------------------------------
% Number of D-H points: (256)^2 points
%   with_h   err = max|U32 - u_exact|
%      1     1.321165399303936e-14
%      2     1.321165399303936e-14
%      3     1.321165399303936e-14
%      4     1.321165399303936e-14 
%   err_noh = max|U32 - u_exact| =  1.321165399303936e-14
% ----------------------------------------------------------
%   with_h  err=max |U64 - u_exact|
%      1    4.329869796038111e-15 
%      2    4.329869796038111e-15 
%      3    4.329869796038111e-15
%      4    4.329869796038111e-15 
%   err_noh = max|U64 - u_exact| =  4.329869796038111e-15       
