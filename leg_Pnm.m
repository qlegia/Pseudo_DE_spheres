% program leg_Pnm
% Sets up the P matrix containing normalised asso. Leg. fns.  P_n^|m|(theta)
% Here theta are the Gauss-quadrature points in [-1,1]
% Note that the we have normalised P_n^|m| in P, rather than just P_n^m

% Uses the matlab files : fac.m (written by us) and legendre.m

function [P] = leg_Pnm(sph_N,theta)
P = [];
for n = 1:sph_N-1
% First we get P_n^m in temp
temp = legendre(n,theta);

% Make temp to contain P_n^|m|
temp = [temp([n+1:-1:2],:);temp];

% Get the normalising factor c_n^|m|
c=(fac(n*ones(2*n+1,1)- [[n:-1:1]';[0:n]']))./(fac(n*ones(2*n+1,1)+[[n:-1:1]';[0:n]']));
c= sqrt(c*(2*n+1)/(4*pi))*ones(1,length(theta));

% Normalise elements in temp
temp = temp.*c;
P = [P;temp];
end
c = [];
temp = [];

