function S = Srand(d)
% create d random points on unit sphere in R^3
% S is 2 by d matrix of spherical coordinates with
% columns corresponding to different points.
% First row is angle theta, second row is angle phi

S(1,1:d) = 2*pi*rand(1,d);      % Theta in [0, 2*pi)
S(2,1:d) = pi*rand(1,d);	    % Phi in [0, pi]
S(1:2,1) = [0; 0];              % First point has phi = 0
S(2,2) = 0;                     % Second point has theta = 0


