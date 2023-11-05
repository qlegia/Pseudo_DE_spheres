function y=legendre_vals(L,z,varargin)
% LEGENDRE_VALS takes as input L and a vector z, 
% where -1<=z(k)<=1, and it returns all
% Schmidt normalized Legendre functions for orders
% l=0 to l=L evaluated at z. The output is a multi-
% dimensional array, y(m+1,k,l+1)=S^m_l(z(k)). NOTE:
% the vector z is turned into a row. This may require
% reshaping when the results are used in other 
% computations. Also, we are calculating associated 
% Legendre polynomials, not spherical harmonics. 
% Thus, for a given l we calculate l+1 values, 
% starting with m=0 and going up to m=l.
%
% An optional argument allows specifiying ARRAY or 
% ROW form for output. Default is ARRAY.

% Put z in the form of a row vector.
x=z(:)';

% Initialize a matrix to be used in computing the weights
% for the recursion formula used below.
w=((1:L)').^2*ones(1,L-1) - ones(L,1)*(0:L-2).^2;

% We only need the lower triangular part. Even so,
% it's convenient to make the upper triangular part 
% all ones, because we will have to do element by element
% division later.
w=tril(w)+triu(ones(L,L-1),1);

% Form the squares of the entries in the weight matrix for the l+1 term.
w1=1./w(1:L-1,:);

% Form the squares of the entries in the weight matrix for the l+2 term.
w2=w(2:L,:).*w1;

% Calculate square roots.
w1=sqrt(w1); w2=sqrt(w2);

% Initialize the 3D array. 
y=zeros(L+1,length(x),L+1);

% Reinializatize matrices (patches bug) and compute low ``L'' values.
ll=max(mod(L,15),1);
for kk=L:-15:ll,
  y(1:kk+1,:,kk+1)=legendre(kk,x,'sch');
  y(1:kk,:,kk)=legendre(kk-1,x,'sch');

  for l=(kk-2):-1:max(kk-15,0)
    m=0:l;
    y(m+1,:,l+1)=diag((2*l+3)*w1(l+1,1+m))*y(m+1,:,l+2)*diag(x)...
	- diag(w2(l+1,1+m))*y(m+1,:,l+3);
  end
end

if nargin>2,
    vn=varargin{1};
    switch lower(vn)
    case 'array'
        y;
    case 'row'
        y_row=[];
        for l=0:L,
            y_row=cat(1,y_row,y(1:l+1,:,l+1));
        end
        y=y_row;
    otherwise
        error('Not row or array')
    end
else
    y;
end
