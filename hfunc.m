function y=hfunc(q,x)
%Call y=hfunc(q,x)
%H. N. Mhaskar December 25, 2003
%evaluates a spline function of order q (polynomials of degree q-1, support on [0,q]
%and smoothness q-2) that is 1 on [0,1/2] and zero outside [0,1].

y=zeros(size(x));
twoq=2*q;
scaledx=twoq.*x;
for j=-q:q
   y=y+bspline(q,scaledx-j);
end
