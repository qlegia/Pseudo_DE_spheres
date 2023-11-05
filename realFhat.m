function FH = realFhat(k,Lmax)
% 15-Sep-2005: first draft
% compute real Fourier coefficients hat{f}(ell,m) for all ell<Lmax
% based on Driscoll-Healy quadrature rule
%------------------------------------------------------------------------------
% Q.T. Le Gia, UNSW, Sydney, Australia. 
% Assume n = 2^k = b^2, Lmax < b/2
%-------------------------------------------------------------------------------
  b = 2^k;
  if (Lmax>=b/2)
      sprintf('We require Lmax < %d',b/2)
      return
  end
  tic
  a = zeros(1,b);
  for j=0:b-1    
     v_ell = 2*[0:b/2-1]+1;    
     a(j+1) = (2*sqrt(2)/b)*sin(pi*j/b)*sum(sin(v_ell*pi*j/b)./v_ell);
  end
  sum(a)
  theta = (pi/b)*[0:1:b-1]';
  phi = (2*pi/b)*[0:1:b-1];
  sintheta = sin(theta);
  X = sintheta*cos(phi);
  Y = sintheta*sin(phi);
  Z = cos(theta)*ones(1,b);
  f = zeros(b,b);
  fh_0 = zeros(1,b);
  fh_plus = zeros(b,b);
  fh_minus = zeros(b,b);
  Cphi = cos((1:b)'*phi); % matrix of cos(m phi_j)  j=0..b-1
  Sphi = sin((1:b)'*phi); % matrix of sin(m phi_j)
  % compute the values of f at sample points
  for k=0:b-1
    % f(k+1,:) = f(theta_k, phi_j) j=0..b-1
    in1 = find(X(k+1,:)> 0.9);
    f(k+1,in1) = f(k+1,in1) + (X(k+1,in1)-0.9).^(3/4);
    in2 = find(Z(k+1,:)>0.9);
    f(k+1,in2) = f(k+1,in2) + (Z(k+1,in2)-0.9).^(3/4);
  end
  % finish compute the values of f at sample points
  for k=0:b-1
    fh_0(k+1) = sum(f(k+1,:)); % for m=0
    fh_plus(:,k+1)  = Cphi*f(k+1,:)'; % for m>0
    fh_minus(:,k+1) = Sphi*f(k+1,:)'; % for m<0 
  end
  Xz = cos(theta);
  % ell=0, m=0 is treated separately
  % legendre(0,x,'sch') = 1 for any x
  FH(1) = sum(a.*fh_0)/sqrt(4*pi);
  for ell=1:Lmax
    % all values of associated P^m_ell(cos(theta_j)), |m|<=ell
    Pell = legendre(ell,Xz,'sch'); 
    r_indx = 1+ell*(ell+1);
    % We want for m>0 
    %  fh_plus(m+1) = sum_{j=0}^{b-1} f(.,phi_j) cos(m phi_j)
    % and for m<0
    %  fh_minus(-m+1) = sum_{j=0}^{b-1} f(.,phi_j) sin(m phi_j) 
    %
    % treat m=0 separately
    FH(r_indx) = sum(a.*Pell(1,:).*fh_0);
    for m=1:ell
      FH(r_indx+m) = sum(a.*Pell(m+1,:).*fh_plus(m,:));
      FH(r_indx-m) = sum(a.*Pell(m+1,:).*fh_minus(m,:));
    end
    FH(r_indx-m:r_indx+m) = sqrt((2*ell+1)/(4*pi))*FH(r_indx-m:r_indx+m);
    ell
    toc
  end
  FH = (FH*sqrt(2*pi)/b)*sqrt(4*pi);
