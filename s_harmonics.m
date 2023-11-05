function Y=s_harmonics(L,theta,phi,type)
%S_HARMONICS computes spherical harmonics of the form (cf. A. R. Edmonds,
%Angular Momentum in Quantum Mechanics, 2nd ed., Princeton U. Press, 1960)
%
%   Y_{\ell,m} = ((2\ell+1)/4pi)S^m_\ell(\cos(\theta)e^{i m \phi}c_m
%
%where c_0=1 and c_m=2^{-1/2}(-1)^m if m>0 and c_m=2^{-1/2} if m<0. 
%The S^m_\ell are the Schmidt normalized associated Legendre functions. 
%The inputs are L, the highest order of the spherical harmonic desired, 
%the co-latitude theta, and the azimuthal angle phi. The output Y is an
%(L+1)^2 x N array. The row index labels the spherical harmonic
%Y_{\ell,m} as r_index=\ell(\ell+1)+m. The column labels points
%given by (theta,phi) spherical coordinates.
%
%If real-valued spherical harmonics are desired, then the variable 'type'
%is set equal to 'real'. In that case, m<0 corresponds to sin(m\phi) and
%m>0 corresponds to cos(m\phi). The variable 'type' can also be set equal 
%to 'complex', which returns the default defined above.
%
%
%Francis J. Narcowich, 17 December 2001.

%Turn theta and phi into rows.
theta=theta(:)';
phi=phi(:)';

%Compute the Schmidt normalized associated Legendre functions of
%order L and less. The 'array' argument means the output s is an
%(L+1) x N x (L+1) dimensional array, s, for which
%s(m+1,k,\ell+1)=S^m_\ell(cos(theta(k)).
s=legendre_vals(L,cos(theta),'array');

if nargin==3,
  type='complex';
end

switch type
  case {'complex'}
    %The matrix phase contains 2^{-1/2}e^{i m \phi} for m>0.
    phase=2^(-1/2)*exp(i*(1:L)'*phi);
    plus_minus=(-1).^(1:L);
    %Initialize the first row of Y.
    Y(1,:)=s(1,:,1)/sqrt(4*pi); %sqrt(4*pi) normalizes l=0 case.
    %Compute the remaining rows.
    for k=1:L,
      y=s(2:(k+1),:,k+1).*phase(1:k,:); %m=1:\ell
      y=[flipud(conj(y));s(1,:,k+1);diag(plus_minus(1:k))*y];%m=-\ell:ell 
      y=sqrt((2*k+1)/(4*pi))*y; %Gives correct normalization for Y_{\ell,m}
      Y=[Y;y]; % Add rows to previous matrix.
    end
  case {'real'}
    C=cos((1:L)'*phi); S=sin((1:L)'*phi);
    Y(1,:)=s(1,:,1)/sqrt(4*pi); %sqrt(4*pi) normalizes l=0 case.
    for k=1:L,
      y_C=s(2:(k+1),:,k+1).*C(1:k,:); %m=1:\ell
      y_S=flipud(s(2:(k+1),:,k+1).*S(1:k,:));%m=-\ell:-1
      y=[y_S;s(1,:,k+1);y_C]*sqrt((2*k+1)/(4*pi)); %Normalize
      Y=[Y;y]; %Add rows to previous matrix.
    end
end
