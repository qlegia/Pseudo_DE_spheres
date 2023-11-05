function X = s2x(S)
% Convert the spherical representation S of a fundamental system
% on the sphere in R^3 into cartesiaon coordinates X
% S is a row vector containing the angles phi and theta
% Because of the rotational invariance of the fundamental system
% one point is always (x,y,z) = (0, 0, 1) <=> phi = 0 (and theta = 0)
% a second point always has theta = 0.
% Thus there are d-1 variables phi and d-2 variables theta
% X is a 3 by d matrix containing the d cartesian points on the sphere
lS = length(S);
d = (lS+3) / 2;
phi = [0, S(1:d-1)];
theta = [0, 0, S(d:lS)];
X = zeros(3,d);
X(1,:) = sin(phi) .* cos(theta);
X(2,:) = sin(phi) .* sin(theta);
X(3,:) = cos(phi);