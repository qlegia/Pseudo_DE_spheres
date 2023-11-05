function FH = FFT_Fhat(k,Lmax)
% 5-Sep-2005: first draft
% 7-Sep-2005: added line 46: if ell>0
% 8-Sep-2005: changed line 38 to legendre(ell,Xz,'sch') to get rid of q_plus
%             and q_minus and the use of fac.m
% 9-Sep-2005: modified line 56,57 to take care factor 1/sqrt(2)
% Based on Theorem 8 in 
%    Driscoll-Healy's Adv. in Math. 15, pp. 202-250 (1994) 
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
  theta = (pi/b)*[0:1:b-1]';
  phi = (2*pi/b)*[0:1:b-1];
  sintheta = sin(theta);
  X = sintheta*cos(phi);
  Y = sintheta*sin(phi);
  Z = cos(theta)*ones(1,b);
  f = zeros(b,b);
  for k=0:b-1
    % f(k+1,:) = f(theta_k, phi_j) j=0..b-1
    Pts = [X(k+1,:); Y(k+1,:); Z(k+1,:)]';
    % change the following line for different function f
    f(k+1,:) = (non_smooth_f(Pts))';
    fh_0(k+1) = sum(f(k+1,:));            % for m=0
    fh_plus(k+1,:) = fft(f(k+1,:),b);     % for m>0
    fh_minus(k+1,:) = b*ifft(f(k+1,:),b); % for m<0 
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
    %   fh_plus(m+1) = sum_{j=0}^{n-1} f(.,phi_j)exp(-i*2*pi*m*j/n)
    % and for m<0
    %   fh_minus(-m+1) = sum_{j=0}^{n-1} f(.,phi_j)exp(i*2*pi*(-m)*j/n) 
    %
    % treat m=0 separately
    FH(r_indx) = sum(a.*Pell(1,:).*fh_0);
    plus_minus = (-1).^(1:ell);
    for m=1:ell
      FH(r_indx+m) = (1/sqrt(2))*plus_minus(m)*sum(a.*Pell(m+1,:).*(fh_plus(:,m+1))');
      FH(r_indx-m) = (1/sqrt(2))*sum(a.*Pell(m+1,:).*(fh_minus(:,m+1))');
    end
    FH(r_indx-m:r_indx+m) = sqrt((2*ell+1)/(4*pi))*FH(r_indx-m:r_indx+m);
    ell
    toc
  end
  FH = FH*sqrt(2*pi)*2*sqrt(pi)/b;
