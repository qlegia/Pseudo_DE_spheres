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
            % integrate poly upto deg 64
X = X_W(:,1:3);
% f = (X(:,1).*cos(X(:,3)) - X(:,3).*sin(X(:,3)) + 0.5*cos(X(:,3))).*exp(X(:,1));
f = non_smooth_f(X);
FH = Fhat(X_W,f',63);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_q = 7;
err_1 = zeros(1,max_q);
err_2 = zeros(1,max_q);
err_3 = zeros(1,max_q);
err_4 = zeros(1,max_q);
%[err_1,errnoh1] = cap_approx_exact(FH,30,XR1,non_smooth_f(XR1'))
%toc
[err_2,errnoh2] = cap_approx_exact(FH,60,XR2,non_smooth_f(XR2'))
toc
%[err_3,errnoh3] = cap_approx_exact(FH,30,XR3,non_smooth_f(XR3'))
%toc
%[err_4,errnoh4] = cap_approx_exact(FH,30,XR4,non_smooth_f(XR4'))
%sprintf('%0.10e %0.10e %0.10e %0.10e \n',errnoh1,errnoh2,errnoh3,errnoh4)
%for q=1:max_q
%  sprintf('%0.10e %0.10e %0.10e %0.10e \n',err_1(q),err_2(q),err_3(q),err_4(q))
%end  
toc