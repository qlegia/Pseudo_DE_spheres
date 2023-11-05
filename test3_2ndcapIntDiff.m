% 2-Sep-2005.
% compute good solutions of degree 100 for integral and differential equations
% and store them as reference solution.
clear all;
tic
load Fhat128; %used Driscoll-Healy rule DH128 
[goodsol,goodsoln_noh] = cap_int_pwh(goodFH,120);
save SolnsInt120 goodsol goodsoln_noh;
toc
tic
load Fhat128;
[goodsol,goodsoln_noh] = cap_diff_pwh(goodFH,120);
save SolnsDiff120 goodsol goodsoln_noh;
toc
