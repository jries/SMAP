function [sx,sy]=getsxfromzfitpar(z,zfit,midp)
z=z/1000;
midp=midp/1000;
%zpar=[sigma0x,Ax,Ay,Bx,By,gamma,d,sigma0y)
fitp=[zfit([7 1 8 2 3 4 5 6]) midp];
[sx,sy]=sbothfromsigma(fitp,z,true);
end

function s=sigmafromz(par,z,B0)
% global gamma
% parx= [d sx0 Ax Bx g mp]
s0=par(2);d=par(1);A=par(3);B=par(4)*B0;g=par(5);mp=par(6);

s=s0*sqrt(1+(z-g+mp).^2/d^2+A*(z-g+mp).^3/d^3+B*(z-g+mp).^4/d^4);
end

function [sx,sy]=sbothfromsigma(par,z,B0)
% parx= [d sx0 sy0 Ax Ay Bx By g mp]
px=par([1 2 4 6 8 9]);
py=par([1 3 5 7 8 9]);
 py(5)=-py(5);
% zh=z(:,1);
sx=sigmafromz(px,z,B0);sy= sigmafromz(py,z,B0);
end