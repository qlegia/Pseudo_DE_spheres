function int_f = mon_carlo(N,f)
% N random points on the unit sphere
% f = vector [1..N] containing [f(x_j): j=1..N]
% Perform Monte-Carlo integration
% Reference: http://mathworld.wolfram.com/MonteCarloIntegration.html
% 1) evaluate <f> 
% 2) evaluate <f^2>
% 3) int(f) by Monte-Carlo method
mf    = sum(f)/N;
msqrf = sum(f.^2)/N;
int_f = 4*pi*mf + 4*pi*sqrt( (msqrf - mf^2)/N);
