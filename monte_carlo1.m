function int_f = monte_carlo1(N,Y,f)
% N random points on the unit sphere
% f is a row vector [1..N] containing [f(x_j): j=1..N]
% Y is a matrix with mth-row = [Y_m(x_1) ... Y_m(x_N)]
% Output is a vector with mth_element = Q( f Y_m)
%  where Q is the Monte-Carlo integration
% -----------------------------------------------------------------
% Perform Monte-Carlo integration
% Reference: http://mathworld.wolfram.com/MonteCarloIntegration.html
% 1) evaluate <f> 
% 2) evaluate <f^2>
% 3) int(f) by Monte-Carlo method
nr = size(Y,1); % number of rows of Y
for k=1:nr
  Yf(k,:) = Y(k,:).*f;
end
Yf = Yf';
mf    = sum(Yf)/N; % sum over all points x_j
msqrf = sum(Yf.^2)/N;
int_f = 4*pi*mf + 4*pi*sqrt( (msqrf - mf.^2)/N);
