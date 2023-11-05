function y=powhfunc(m,x)
%Written by Dr. H. N. Mhaskar, August 30, 2005
%y=powhfunc(m,x)
%evaluates the function $h(x)=1-A\int_{1/2}^x (t-1/2)^m(1-t)^m dt$,
%where A is chosen so that $h(1)=0$.
%To be used only for $m\le 9$, otherwise the matlab factorial is not accurate.
%Also, x is assumed to be in the range [0,1].

factor=(factorial(2*m+1)/factorial(m)).*((2.*x-1).^(m+1));
sum=zeros(size(x));
for k=0:m
   sum=sum+((1-2.*x).^k)./(factorial(k)*factorial(m-k)*(m+k+1));
end
indless=find(x<=1/2);
indhigh=find(x>1/2);
y(indless)=1;
y(indhigh)=1-factor(indhigh).*sum(indhigh);