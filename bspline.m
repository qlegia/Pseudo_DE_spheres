function y=bspline(q,x)
%call bspline(q,x)
%Compute a cardinal B-spline of order $q$ (piecewise polynomials of degree
% $q-1$, support $[0,q]$) using the recurrence in C. de Boor: "A practical guide
%to splines", Springer Verlag, New York, 1978, p. 131.

y=zeros(size(x));
if (q==1)
   I= find((x>0)& (x <=1));
   y(I)=1;
  else 
   y = (x.*bspline(q-1,x)+(q-x).*bspline(q-1,x-1))/(q-1);
end  
