function zfit=getzfitpar(sx,sy,znm,zrange,midpoint,B0,ax)

% mpf=midpreal-midp;
% parx= [d sx0 sy0 Ax Ay Bx By g mp]
startp=[    0.3    1.0    1.0000  0   0        0         0  0.307   -midpoint/1000];
% startp=[    0.3    1.0    1.0000  0   0        0         0  0.307   -mpf];

startp(2)=myquantile(sx,0.01);startp(3)=myquantile(sy,0.01);
% startp(3)=1.5;
ind=znm>zrange(1)&znm<zrange(2);
sx=sx(ind);
sy=sy(ind);
znm=znm(ind);
z=znm/1000;

% B0=false;
fitp=lsqnonlin(@sbothfromsigmaerr,startp,[],[],[],[z z],[sx sy],0);
if B0
fitp=lsqnonlin(@sbothfromsigmaerr,fitp,[],[],[],[z z],[sx sy],true);
end

zt=min(z):0.01:max(z);
if nargin>6
    
    
    sxfs=sigmafromz(startp([1 2 4 6 8 9]),zt,B0);
sxf=sigmafromz(fitp([1 2 4 6 8 9]),zt,B0);
startpy=startp([1 3 5 7 8 9]);
startpy(5)=-startpy(5);
syfs=sigmafromz(startpy,zt,B0);

plot(ax,z,sx,'r.')
hold on
plot(ax,z,sy,'r.')
plot(ax,zt,sxf,'k')
% plot(ax,zt,sxfs,'c',zt,syfs,'g')
fpy=fitp([1 3 5 7 8 9]);
fpy(5)=-fpy(5);
syf=sigmafromz(fpy,zt,B0);
plot(ax,zt,syf,'k')
% hold off
axis tight
ylabel(ax,'PSFx,PSFy')
end
%zpar=[sigma0x,Ax,Ay,Bx,By,gamma,d,sigma0y)
zfit=real(fitp([2 4 5 6 7 8 1 3]));
end

function s=sigmafromz(par,z,B0)
% global gamma
% parx= [d sx0 Ax Bx g mp]
s0=par(2);d=par(1);A=par(3);B=par(4)*B0;g=par(5);mp=par(6);

s=s0*sqrt(1+(z-g+mp).^2/d^2+A*(z-g+mp).^3/d^3+B*(z-g+mp).^4/d^4);
end

function s=sbothfromsigma(par,z,B0)
% parx= [d sx0 sy0 Ax Ay Bx By g mp]
px=par([1 2 4 6 8 9]);
py=par([1 3 5 7 8 9]);
 py(5)=-py(5);
zh=z(:,1);
s=[sigmafromz(px,zh,B0) sigmafromz(py,zh,B0)];
end

function err=sbothfromsigmaerr(par,z,sx,B0)
sf=sbothfromsigma(par,z,B0);
err=sf-sx;

% stderr=std(err);
err=err./sqrt(abs(err));
end