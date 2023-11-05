% To get a good solution, then the good solution is saved to "goodU32.mat"
% -------------------------------------------------------------------------
% Oct-04-2004: use new_spectral
% -------------------------------------------------------------------------
clear all;
load DH_64; % change this command to load different data sets
X = X_W(:,1:3);
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
[FH,u32] = new_spectral(X_W,f',50,32,1);
save goodU32 u32;
