% Spectral method using Monte-Carlo integration
% Q. T. Le Gia, UNSW, Australia, 22-Aug-2005.
% Sep-19-2005: fixed bug X=X' and rerun.
% --------------------------------------------------
clear all;
d = 65536;
S = srand(d); % generates 1000 random points on S^2
X = s2x([S(1,2:d),S(2,3:d)]);   % convert to Cartesian points
X = X'; %make sure X has 3 columns
f = (X(:,1).*cos(X(:,3)) - X(:,3).*sin(X(:,3)) + 0.5*cos(X(:,3))).*exp(X(:,1));
[FH,sol] = MCspectral(X,f',64,4);
% now we have to compute the error between the exact solution and 
% the approximation
num = 50;
[Sx,Sy,Sz] = sphere(num);
err_max = zeros(num+1,num+1);
for k=1:num+1       
       err_max(k,:) = abs(sol(k,:) - exp(Sx(k,:)).*cos(Sz(k,:)));               
end    
max(max(err_max)) 
norm(err_max)/(num+1)
% Results after fixed X=X' Sep-19
%   65536 points
%    n = Lsol = 64
%         1st run           2nd run           3rd run
%      m  max_err l2_err    max_err l2_err    max_err  l2_err
%      0  1.6722  0.7065    1.7504  0.7135    1.7315   0.7215
%      1  2.1575  1.1313    2.1590  1.1334    2.1148   1.1288
%      2  2.0400  1.1240    1.7408  0.7195    1.7470   0.7235
%      3  2.0675  1.1268    1.7030  0.7162    1.7222   0.7149    
%      4  2.0498  1.1268    1.7075  0.7277    1.6872   0.7210 
%      
