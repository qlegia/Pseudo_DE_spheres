% program fac
% computes the factorial of each element in a vector "m"

function [val] = fac(m)
if(size(m,2) > 1)
m = m';
end
n = max(m);
s = length(m);
x = ones(s,1)*(1:n);
%[i,j] = find(((x - m*ones(1,n)) > 0))   %double index has trouble, use single
jj = find((x - m*ones(1,n)) > 0);
%x(i,j) = 1;                             %double index has trouble, use single 
x(jj) = 1;
val = (prod(x'))';
