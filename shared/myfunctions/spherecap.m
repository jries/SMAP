function [x,y,z,n]=spherecap(par,dx)
% dx: distance in real coordinates
R=par(1);
angle=par(2);
da=dx/R;

theta=da/2:da:angle+da/2;


numphi=round(2*pi*R/dx);
phi=linspace(0,2*pi,numphi);
phi(end)=[];

[T,P]=meshgrid(theta,phi);

x=R*sin(T).*cos(P);
y=R*sin(T).*sin(P);
z=R*cos(T);
n=sin(T);
