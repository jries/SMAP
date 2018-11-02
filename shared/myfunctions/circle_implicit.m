function error=sphere_implicit(par,x,y,z,wx,wy,wz)
%par(1)=r; par(2:4)=x0,y0,z0
% error2=(x-par(2)).^2+(y-par(3)).^2+(z-par(4)).^2-par(1)^2;

error=sqrt(((x-par(2)).^2)+((y-par(3)).^2))-par(1);