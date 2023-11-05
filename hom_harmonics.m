function Y=hom_harmonics(L,theta,phi,type)
%HOM_HARMONICS computes homogeneous spherical harmonics of the form 
%(cf. A. R. Edmonds, Angular Momentum in Quantum Mechanics, 2nd ed., 
% Princeton U. Press, 1960)
%
%   Y_{\ell,m} = (sqrt(2\ell+1)/4pi))S^m_\ell(\cos(\theta)e^{i m \phi}c_m
%
%where c_0=1 and c_m=2^{-1/2}(-1)^m if m>0 and c_m=2^{-1/2} if m<0. 
%The S^m_\ell are the Schmidt normalized associated Legendre functions. 
%The inputs are L, the order of the spherical harmonic desired, 
%the co-latitude theta, and the azimuthal angle phi. The output Y is an
%(2L+1) x N array. The row index labels the spherical harmonic
% Y_{L,m} for m=-L...L as [Y_{L,-L},...,Y_{L,0},...,Y_{L,L}]. The column 
% labels points given by (theta,phi) spherical coordinates.
%
%If real-valued spherical harmonics are desired, then the variable 'type'
%is set equal to 'real'. In that case, m<0 corresponds to sin(m\phi) and
%m>0 corresponds to cos(m\phi). The variable 'type' can also be set equal 
%to 'complex', which returns the default defined above.

%Turn theta and phi into rows.
theta=theta(:)';
phi=phi(:)';

%Compute the Schmidt normalized associated Legendre functions of order L. 
%s(m+1,k)=S^m_L(cos(theta(k)).
s=legendre(L,cos(theta),'sch');

if nargin==3,
  type='complex';
end

switch type
  case {'complex'}
    %The matrix phase contains 2^{-1/2}e^{i m \phi} for m>0.
    phase=2^(-1/2)*exp(i*(1:L)'*phi);
    plus_minus=(-1).^(1:L);
    %Initialize the first row of Y.
    if (L==0)
       Y(1,:)=s(1,:)/sqrt(4*pi); %sqrt(4*pi) normalizes l=0 case.
    else       
       y=s(2:L+1,:).*phase(1:L,:); %m=1:L
       Y=[flipud(conj(y));s(1,:);diag(plus_minus(1:L))*y];%m=-L:L 
       Y=sqrt((2*L+1)/(4*pi))*Y; %Gives correct normalization for Y_{L,m}
    end
  case {'real'}
    C=cos((1:L)'*phi); S=sin((1:L)'*phi);
    if (L==0)
       Y(1,:)=s(1,:,1)/sqrt(4*pi); %sqrt(4*pi) normalizes l=0 case.
    else 
      y_C=s(2:(L+1),:).*C(1:L,:); %m=1:L
      y_S=flipud(s(2:(L+1),:).*S(1:L,:));%m=-L:-1
      Y=[y_S;s(1,:);y_C]*sqrt((2*L+1)/(4*pi)); %Normalize      
    end
end