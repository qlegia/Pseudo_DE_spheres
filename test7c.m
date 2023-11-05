% In order to use this routine, you have to run "test7_U39" first 
% to get a good solution, then the good solution is saved to "int_goodU59.mat"
% -------------------------------------------------------------------------
% Nov-15-2004: use integral2
% -------------------------------------------------------------------------
clear all;
tic
load md27.0784; % change this command to load different data sets
X = md27(:,1:3);
zthres = 0.9;
xthres = 0.9;
in1 = find(X(:,3) > zthres);
in2 = find(X(:,1) > xthres);
f = zeros(size(X(:,3)));
f1 = f;
f2 = f;
f1(in1) = X(in1,3)-zthres;
f2(in2) = X(in2,1)-xthres;
f = f1+f2;
f=f.^(3/4);
[Sx,Sy,Sz] = sphere(50);
in3 = find(Sz > zthres);
in4 = find(Sx > xthres);
fplot = zeros(size(Sz));
fplot(in3) = Sz(in3) - zthres;
fplot(in4) = fplot(in4) + (Sx(in4) - xthres);
fplot = fplot.^(3/4);
%surf(Sx,Sy,Sz,fplot);
%colorbar('vert');
%shading interp;
[FH,soln] = integral2(md27,f',27,13,9);
load int_goodU22;
err = abs(soln-u22);
% compute local error on a spherical cap C := {x . w>0.9}
% with w = (1/sqrt(2),0,1/sqrt(2))
max_e_cap = cap_err(51,err,[1/sqrt(2);0;1/sqrt(2)],0.9)
max_e_cap2 = cap_err(51,err,[-1/sqrt(2);0;-1/sqrt(2)],0.9)
max_e_cap3 = cap_err(51,err,[-1/sqrt(2);0;-1/sqrt(2)],0.93)
max_e_cap4 = cap_err(51,err,[-1/sqrt(3);-1/sqrt(3);-1/sqrt(3)],0.9)
max_err = max(max(err)) 
%err = err/max_err;
%figure;
%surf(Sx,Sy,Sz,err);
%colorbar('vert');
%shading interp;
toc