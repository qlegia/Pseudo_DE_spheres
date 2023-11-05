% 31-Aug-2005;
clear all;
tic
c2 = [-1/sqrt(2);0;1/-sqrt(2)];          cos_al2 = 0.9;
alpha2 = acos(cos_al2);
XR2 = ranpts_cap(1000,c2,alpha2);
[theta,phi]=x2theta(XR2);
Y0_100 = s_harmonics(100,theta,phi,'real'); %54.25 seconds
save Y0_100_XR2 Y0_100 XR2;
toc
Y101_150 = [];
for ell=101:150
  Y = hom_harmonics(ell,theta,phi,'real');
  Y101_150 = [Y101_150;Y];
end
save Y101_150_XR2 Y101_150 XR2; % 111.72 seconds
% load Driscoll-Healy quad rule
% load DH_64;
% load md191.36864;
% X = md191(:,1:3);
% f = non_smooth_f(X);
% goodFH = Fhat(md191,f',180);
% save Fhat180 goodFH; 
% Fhat128 uses DH_128 takes 1882.82 seconds 
% Fhat64  uses DH_64  takes  109.63 seconds
% Fhat180 uses md191  takes 1932.83 seconds 
toc
