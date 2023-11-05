% test_cap_spectral
% 30-Aug-2005
% Q.T. Le Gia, UNSW, Sydney, Australia.
clear all;
tic
% we want solution on a spherical cap centered at ctradius alpha
Nr = 1000; %number of random points for ssolutions 

c1 = [1/sqrt(2);0;1/sqrt(2)];            cos_al1 = 0.9;
c2 = [-1/sqrt(2);0;1/-sqrt(2)];          cos_al2 = 0.9;
c3 = c2;                                 cos_al3 = 0.93;
c4 = [-1/sqrt(3);-1/sqrt(3);-1/sqrt(3)]; cos_al4 = 0.9;

alpha1 = acos(cos_al1);
alpha2 = acos(cos_al2);
alpha3 = acos(cos_al3);
alpha4 = acos(cos_al4);

XR1 = ranpts_cap(Nr,c1,alpha1);
XR2 = ranpts_cap(Nr,c2,alpha2);
XR3 = ranpts_cap(Nr,c3,alpha3);
XR4 = ranpts_cap(Nr,c4,alpha4);
% compute Fourier coefficients for the good solution
load DH_64; % Driscoll-Healy quadrature X_W of 2^{14} points
X = X_W(:,1:3);
% f = (X(:,1).*cos(X(:,3)) - X(:,3).*sin(X(:,3)) + 0.5*cos(X(:,3))).*exp(X(:,1));
f = non_smooth_f(X);
goodFH = Fhat(X_W,f',63);
% compute Fourier coeffs for the approximation
clear X_W;
clear X;
clear f;
load DH_32;
X = X_W(:,1:3);
%f = (X(:,1).*cos(X(:,3)) - X(:,3).*sin(X(:,3)) + 0.5*cos(X(:,3))).*exp(X(:,1));
f = non_smooth_f(X);
FH = Fhat(X_W,f',31);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_q = 7;
err_1 = zeros(1,max_q);maxg1 = zeros(1,max_q);maxs1 = zeros(1,max_q);
err_2 = zeros(1,max_q);maxg2 = zeros(1,max_q);maxs2 = zeros(1,max_q);
err_3 = zeros(1,max_q);maxg3 = zeros(1,max_q);maxs3 = zeros(1,max_q);
err_4 = zeros(1,max_q);maxg4 = zeros(1,max_q);maxs4 = zeros(1,max_q);
[err_1,maxg1,maxs1] = cap_approx(goodFH,60,FH,20,XR1)
toc
[err_2,maxg2,maxs2] = cap_approx(goodFH,60,FH,20,XR2)
[err_3,maxg3,maxs3] = cap_approx(goodFH,60,FH,20,XR3)
[err_4,maxg4,maxs4] = cap_approx(goodFH,60,FH,20,XR4)
for q=1:max_q
  sprintf('%0.10e %0.10e %0.10e %0.10e \n',err_1(q),err_2(q),err_3(q),err_4(q))
end  
toc
