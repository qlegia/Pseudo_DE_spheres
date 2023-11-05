function pts=randpts(no_points,alpha, beta)
%written by H. N. Mhaskar, Dec. 11, 2002
%call pts=randpts(no_points,alpha, beta)
% alpha>beta; to find on cap of radius alpha, use beta=0
%calculate no_points randomly distributed points on the spherical cap
%\{\cos\alpha \le z\le \cos\beta\}.
% returns a no_points\times 3 matrix, giving the cartesian coordinates.

y=zeros(no_points,2);
randvects=rand(no_points,2);
%matlab recommends the following restart of the rand function
rand('state',sum(100*clock));
start=1;
n=floor(sqrt(no_points));
factor=2*n*n*sin((alpha-beta)/(2*n))/(cos(beta)-cos(alpha));
for k=1:n
   %find number of points in the k-th band of phi
   bandpts=floor(factor*sin(beta+(alpha-beta)*(2*k-1)/(2*n)));
   if (bandpts<=0) bandpts=1; end
   %distribute bandpts points in the theta direction.
   y(start:bandpts+start-1,2)=randvects(start:bandpts+start-1,2).*((alpha-beta)/n)...
      +beta+(alpha-beta)*(k-1)/n;
   y(start:bandpts+start-1,1)=randvects(start:bandpts+start-1,1).*(2*pi)-pi;
   start=start+bandpts;
end
if (start<no_points)
   y(start:no_points,2)=randvects(start:no_points,2).*(alpha-beta)+beta;
   y(start:no_points,1)=randvects(start:no_points,1).*(2*pi)-pi;
end


y(:,2)=pi/2-y(:,2);%to adjust for the weird convention of matlab
%Show the points graphically.
[xx,yy,zz]=sph2cart(y(:,1),y(:,2),ones(no_points,1));
pts=[xx,yy,zz];
return

minz=min(z)
maxz=max(z)
[XI, YI,ZI]=sphere(n);
hold
mesh(XI,YI,ZI)
plot3(x,y,z, '.')
hold off
