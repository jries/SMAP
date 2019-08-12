function simulationerror(locgt,locfit,psf,searchradius)

pixelsize=100
if ~isfield(locgt,'x'), locgt.x=locgt.xnm; end
if ~isfield(locfit,'x'), locfit.x=locfit.xnm; end
if ~isfield(locgt,'y'), locgt.x=locgt.ynm; end
if ~isfield(locfit,'y'), locfit.x=locfit.ynm; end
if ~isfield(locgt,'z'), locgt.x=locgt.znm; end
if ~isfield(locfit,'z'), locfit.x=locfit.znm; end

if nargin<4 || isempty(searchradius)
    searchradius=200;
end

[iAa,iBa,nA,nB,nseen]=matchlocsall(locgt,locfit,0,0,searchradius(1));
totallocs=length(locfit.x);
falsepositives=length(nB);
falsenegatives=length(nA);
matched=length(iAa);

dz=locgt.z(iAa)-locfit.z(iBa);
indinz=abs(dz)<searchradius(end);
dz=dz(indinz);
dx=locgt.x(iAa(indinz))-locfit.x(iBa(indinz));
dy=locgt.y(iAa(indinz))-locfit.y(iBa(indinz));


if nargin<3 || isempty(psf)
    [lp,errphot]=Mortensen(locgt.phot(iAa(indinz)),locgt.bg(iAa(indinz)),150, pixelsize,0);
%     crlb=1000./sqrt(locgt.phot(iAa(indinz))); %
    dxr=dx./lp; %xeems to be closer to 1
    dyr=dy./lp;
    dzr=dz./lp/3; 
else    
    crlb=psf.crlb(locgt.phot(iAa(indinz)),locgt.bg(iAa(indinz)),locgt.znm(iAa(indinz)));
    dxr=dx./sqrt(crlb(:,2))/pixelsize; %xeems to be closer to 1
    dyr=dy./sqrt(crlb(:,1))/pixelsize;
    dzr=dz./sqrt(crlb(:,5));
end

figure(88);
subplot(3,3,1)
fithistr(dx,1)
xlabel('dx')
subplot(3,3,2)
fithistr(dy,1)
xlabel('dy')
subplot(3,3,3)
fithistr(dz,5)
xlabel('dz')

subplot(3,3,4)
fithistr(dxr,0.2)
xlabel('dx/sqrt(CRLBx)')
subplot(3,3,5)
fithistr(dyr,0.2)
xlabel('dy/sqrt(CRLBy)')
subplot(3,3,6)
fithistr(dzr,0.2)
xlabel('dz/sqrt(CRLBz)')

subplot(3,3,7)
hold off
dscatter(locgt.z(iAa),locfit.z(iBa))
hold on
plot([min(locgt.z(iAa)),max(locgt.z(iAa))],[min(locgt.z(iAa)),max(locgt.z(iAa))],'m')
% ,'.',locgt.z(iAa),locgt.z(iAa),'k')
xlabel('z_gt')
ylabel('z_fit')

subplot(3,3,8);
ff='%2.0f';
title(['fpos: ' num2str(falsepositives/totallocs*100,ff), '%, fneg: ' num2str(falsenegatives/totallocs*100,ff) '%'])
% lgt=copyfields([],locgt,{'x','y','znm','phot','bg'})
%matchlocs in defined region
% false pos, false neg
% matched: histogram dx/crlbx dy/crlby dz/crlbz
% std of these quantities vs phot, vs z
subplot(3,3,9);plot(locfit.x,locfit.y,'.',locgt.x,locgt.y,'.')

end

function fithistr(de,dn)
ff='%1.1f';
qq=quantile(de,[0.002 0.998]);
n=floor(qq(1)):dn:ceil(qq(end));
hold off
histogram(de,n);
hn=histcounts(de,n);
nf=n(1:end-1)+(n(2)-n(1))/2;
fitp=fit(nf',hn','gauss1');
ss=fitp.c1/sqrt(2);

% try fitting with second gauss
mp=find(nf>fitp.b1,1,'first');
dn2=ceil(2*ss/dn);
n2=(mp-dn2:mp+dn2)';

% fitp2=fit(nf(n2)',hn(n2)','gauss2','StartPoint',[fitp.a1,fitp.b1,fitp.c1,fitp.a1/10,fitp.b1,fitp.c1*10]);
fitp2=fit(nf(n2)',hn(n2)','gauss1','StartPoint',[fitp.a1,fitp.b1,fitp.c1]);
ss2=fitp2.c1/sqrt(2);
hold on
plot(nf,fitp(nf),'g')
plot(nf(n2),fitp2(nf(n2)),'r')
ingauss=fitp2.a1*sqrt(pi)*fitp2.c1/length(de)/dn;
% de=de(abs(de)<3);
title([num2str(mean(de),2) '±' num2str(std(de),ff) ', fit: ' num2str(fitp2.b1,ff) '±' num2str(ss2,ff) ', in Gauss ' num2str(ingauss*100,'%2.0f') '%'])
xlim([-5*ss 5*ss])
end


function [lp,errphot]=Mortensen(N,Bg,PSF, pixel,cmosn)
b=sqrt(Bg+cmosn^2);

PSFa=sqrt(PSF^2+pixel^2/12);
v=PSFa^2./N.*(16/9+8*pi*PSFa^2*b.^2./N/pixel^2);
lp=sqrt(v);
errphot=sqrt(PSFa^2./N.*(8*pi*PSFa^2.*b.^2./N/pixel^2));
end