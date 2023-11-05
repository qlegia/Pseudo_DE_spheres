function y=bspline_real(q,x)
%call bspline(q,x)
%Compute a cardinal B-spline of order $q$ (piecewise polynomials of degree
% $q-1$, support $[0,q]$) using the recurrence in C. de Boor: "A practical guide
%to splines", Springer Verlag, New York, 1978, p. 131.
%-------------------------------------------------------------------------
y=zeros(size(x));
lenX = length(x);
step = abs(x(2)-x(1)); % assume equal steps
x0 = min(x);
for j=1:lenX
  ty = bspline_int(q,[x0,x0+1]);
  y(j) = ty(1);
  x0 = x0 + step;
end