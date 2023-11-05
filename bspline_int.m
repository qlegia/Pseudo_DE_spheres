function y=bspline_int(q,x)
%call bspline(q,x)
%Compute a cardinal B-spline of order $q$ (piecewise polynomials of degree
% $q-1$, support $[0,q]$) using the recurrence in C. de Boor: "A practical guide
%to splines", Springer Verlag, New York, 1978, p. 131.
%-------------------------------------------------------------------------
% we assume that x(j+1)-x(j) = 1 for all j
y=zeros(size(x));
if (q==1)
   I= find((x>0)& (x <=1));
   y(I)=1;
  else 
   y = (x.*bspline_int(q-1,x)+(q-x).*bspline_int(q-1,x-1))/(q-1);
end  
