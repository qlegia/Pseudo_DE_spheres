d = 360000;
S = srand(d); % generates 1000 random points on S^2
X = s2x([S(1,2:d),S(2,3:d)]);   % convert to Cartesian points      
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
%---------------------------
f = ones(1,d);
int_f = mon_carlo(d,f)

