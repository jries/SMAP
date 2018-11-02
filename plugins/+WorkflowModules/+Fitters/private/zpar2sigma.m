
function [sx,sy]=zpar2sigma(z,zpar)
% zpar=PSFSigma,Ax,Ay,Bx,By,gamma,d,PSFy0);
if length(zpar)==8
    zpar(9)=0;
end
% parx= [d sx0 sy0 Ax Ay Bx By g mp]
px=zpar([7 1 2 4 6 9]);
py=zpar([7 8 3 5 6 9]);
 py(5)=-py(5);

sx=sigmafromz(px,z); 
sy=sigmafromz(py,z);
end

function s=sigmafromz(par,z)
% global gamma
% parx= [d sx0 Ax Bx g mp]
s0=par(2);d=par(1);A=par(3);B=par(4);g=par(5);mp=par(6);
s=s0*sqrt(1+(z-g+mp).^2/d^2+A*(z-g+mp).^3/d^3+B*(z-g+mp).^4/d^4);
end