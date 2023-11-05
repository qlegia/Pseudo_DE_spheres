function y=hfunc2(q,x)
%Call y=hfunc2(q,x)
%evaluates a spline function of order q 
%(polynomials of degree q-1, support on [0,q]
%and smoothness q-2) that is 1 on [0,1/2] and zero outside [0,1].
%-----------------------------------------------------------------
% Q T Le Gia, 12-Sep-2005.
%-----------------------------------------------------------------
if (q==1)
  y = ones(1,length(x));
else
  y=zeros(size(x));
  idx = find(x<=0.5);
  len = length(idx);
  tail = 1:0.1:2;
  len_t = length(tail);
  % using interpolation routine from the Spline toolbox 
  vecx = [x(idx) tail];
  vecy = [ones(1,len) zeros(1,len_t)];
  sp = spapi(q,vecx,vecy);
  y=fnval(sp,x);
end
