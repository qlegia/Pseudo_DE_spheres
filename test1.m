% MATLAB code used to generate Table 1 in mhaskar_report.tex
% Q. T. Le Gia, Oct-6-2004.
% -----------------------------------------------------------------------------------
clear all;
load md19.0400;
load md39.1600;
%load md49.2500; 
%load md59.3600; % load the set of discrete points from Womersley's quadrature
%--------------------------------------------------------------------------
X = md19(:,1:3);
%X = md39(:,1:3);
%X = md49(:,1:3);
%X = md59(:,1:3);
%--------------------------------------------------------------------------
f = (X(:,1).*cos(X(:,3)) - X(:,3).*sin(X(:,3)) + 0.5*cos(X(:,3))).*exp(X(:,1));
[FH,sol] = new_spectral(md19,f',19,12,0);
%[FH,sol] = new_spectral(md39,f',39,22,0);
%[FH,sol] = new_spectral(md49,f',49,25,0);
%[FH,sol] = new_spectral(md59,f',59,22,0);
% u_exact = exp(X(:,1)).*cos(X(:,3));
% now we have to compute the error between the exact solution and the approximation
num = 50;
[Sx,Sy,Sz] = sphere(num);
err_max = zeros(num+1,num+1);
for k=1:num+1       
       err_max(k,:) = abs(sol(k,:) - exp(Sx(k,:)).*cos(Sz(k,:)));               
end    
max(max(err_max))
norm(err_max)/(num+1)
%-------------------------------------------------------------------------
% with the masking function h
%-------------------------------------------------------------------------
[FH,sol] = new_spectral(md19,f',19,12,1);
%[FH,sol] = new_spectral(md39,f',39,22,1);
%[FH,sol] = new_spectral(md49,f',49,25,1);
%[FH,sol] = new_spectral(md59,f',59,22,1);
% u_exact = exp(X(:,1)).*cos(X(:,3));
% now we have to compute the error between the exact solution and the approximation
num = 50;
[Sx,Sy,Sz] = sphere(num);
err_max = zeros(num+1,num+1);
for k=1:num+1       
       err_max(k,:) = abs(sol(k,:) - exp(Sx(k,:)).*cos(Sz(k,:)));               
end    
max(max(err_max))
norm(err_max)/(num+1)

