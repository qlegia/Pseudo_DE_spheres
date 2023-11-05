% MATLAB code used to generate Table 3 in numerical experiments rep
% Q. T. Le Gia, Oct-28-2004.
% -----------------------------------------------------------------------------------
clear all;
%load me19.0400;
load me39.1600;
%load me49.2500; 
%load me59.3600; % load the set of discrete points from Womersley's quadrature
%--------------------------------------------------------------------------
%X = me19(:,1:3);
X = me39(:,1:3);
%X = me49(:,1:3);
%X = me59(:,1:3);
%------------------------------------------------------------------------
f = 2*ones(length(X),1);
%[FH,sol] = integral2(me19,f',19,12,1);
[FH,sol] = integral2(me39,f',39,20,1);
%[FH,sol] = integral2(me49,f',49,25,1);
%[FH,sol] = integral2(me59,f',59,34,1);
% now we have to compute the error between the exact solution and the approximation
num = 50;
[Sx,Sy,Sz] = sphere(num);
err_max = zeros(num+1,num+1);
for k=1:num+1       
       err_max(k,:) = abs(sol(k,:) - 1);               
end    
max(max(err_max))
norm(err_max)/(num+1)
