% 2-Sep-2005.
% compute solutions of degree 50 for integral and differential equations
% and compare with the reference solution of degree 100.
clear all;
tic
load Fhat96;  
load U100_int;
goodsol_noh = zeros(1,1000); % we haven't computed yet
[err1,err_noh1] = cap_int_pwh1(goodsol,goodsol_noh,Fhat96,50);
toc
tic
load U100_diff;
goodsol_noh = zeros(1,1000); % we haven't computed yet
[err2,err_noh2] = cap_diff_pwh1(goodsol,goodsol_noh,Fhat96,50);
toc
% Results as of 2-Sep-2005
% --------------------------------------------------------------------------
% Test 1  . Used Driscoll-Healy DH128 to compute hat{f}_{ell,m}
%         . Used powhfunc.m
%         . f = (x-0.9)^{3/4}_{+} + (z-0.9)^{3/4}_{+} , (x,y,z) on S^2
%         . Errors measure on a cap C centered [-1/sqrt(2);0;-1/sqrt(2)] 
%  a)   Solving integral quation  L u = f
%        err = max| U50 - U100|
%   1.0e-05*
% 1     0.12220176248603   
% 2     0.03457393532649   
% 3     0.01370335475928
% 4     0.00526091292340   
% 5     0.00372259512077   
% 6     0.00220754593720 
% 7     0.00106636778290
% --------------------------------------------------------------------------
%  b)  Solving differential equation L u = f
%       err = max|U50 - U100|
%   1.0e-07 *
% 1    0.43058001877525   
% 2    0.11420361923503   
% 3    0.04492418751237 
% 4    0.01890809965197   
% 5    0.01161643996272   
% 6    0.00834895744553
% 7    0.00424405897578
% --------------------------------------------------------------------------
% Test 2) replace DH64 as quadrature rule
%  Integral equation
%  err = |U50 - U100|
%1.0e-05 *
%
%0.31580727876291
%   0.22246101872359
%   0.20797896718259
%   0.21079018222507
%   0.21590141545326
%   0.21455182245438
%   0.21167856199564
%
%  Differential equation
%  err = |U50 - U100|
% 1.0e-04 *
%
%   0.11629813244143
%   0.11664125220101
%   0.11669417426884
%   0.11666065870058
%   0.11664155351141
%   0.11664334442990
%   0.11665202604698
ans =

   1.0e-05 *

   0.17704179000090
   0.09116663124253
   0.06953930175673
   0.06875178577086
   0.07175164068167
   0.07036389200381
   0.06747642857729

>> err2'

ans =

   1.0e-05 *

   0.34162220945173
   0.34416440441148
   0.34418029869025
   0.34382557835451
   0.34368677818905
   0.34373074136278
   0.34381751184611

