% This script generates a good solution for the Neumann problem
% using R.Womersley's quadrature rules over 3600 points on the sphere.
%----------------------------------------------------------------------
% Q. T. Le Gia - 3-May-2004 - Texas A & M University, College Station.
% Oct-7-2004, use new_spectral with minimized energy points
%----------------------------------------------------------------------
clear all;
load me59.3600; % change this command to load different data sets
X = me59(:,1:3);
zthres = 0.9;
xthres = 0.9;
in1 = find(X(:,3)>zthres);
in2 = find(X(:,1)>xthres);
f = zeros(size(X(:,3)));
f1=f;
f2=f;
f1(in1) = X(in1,3)-zthres;
f2(in2) = X(in2,1)-xthres;
f = f1+f2;
f = f.^(3/4);
[Sx,Sy,Sz] = sphere(50);
in3 = find(Sz>zthres);
in4 = find(Sx>xthres);
fplot = zeros(size(Sz));
fplot(in3) = Sz(in3)-zthres;
fplot(in4) = fplot(in4) + (Sx(in4)-xthres);
fplot = fplot.^(3/4);
surf(Sx,Sy,Sz,fplot);
colorbar('vert');
shading interp;
[FH,u59] = new_spectral(me59,f',59,22,1);
save goodU59 u59;