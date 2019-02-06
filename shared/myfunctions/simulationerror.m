function simulationerror(locgt,locfit,psf)
pixelsize=100
locfit.x=locfit.xnm;locfit.y=locfit.ynm;
maxd=200;%search radius in x,y
[iAa,iBa,nA,nB,nseen]=matchlocsall(locgt,locfit,0,0,maxd);
totallocs=length(locfit.x);
falsepositives=length(nB);
falsenegatives=length(nA);
matched=length(iAa);
dx=locgt.x(iAa)-locfit.x(iBa);
dy=locgt.y(iAa)-locfit.y(iBa);
dz=locgt.znm(iAa)-locfit.znm(iBa);

crlb=psf.crlb(locgt.phot(iAa),locgt.bg(iAa),locgt.znm(iAa));
dxr=dx./sqrt(crlb(:,2))/pixelsize; %xeems to be closer to 1
dyr=dy./sqrt(crlb(:,1))/pixelsize;
dzr=dz./sqrt(crlb(:,5));

figure(88);
subplot(3,3,1)
fithistr(dx,1)
xlabel('dx')
subplot(3,3,2)
fithistr(dy,1)
xlabel('dy')
subplot(3,3,3)
fithistr(dz,1)
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
plot(locgt.znm(iAa),locfit.znm(iBa),'.',locgt.znm(iAa),locgt.znm(iAa),'k')
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
ff='%1.2f';
qq=quantile(de,[0.002 0.998]);
n=floor(qq(1)):dn:ceil(qq(end));
hold off
histogram(de,n);
hn=histcounts(de,n);
nf=n(1:end-1)+(n(2)-n(1))/2;
fitp=fit(nf',hn','gauss1');
ss=fitp.c1/sqrt(2);
hold on
plot(nf,fitp(nf),'r')
ingauss=fitp.a1*sqrt(pi)*fitp.c1/length(de)/dn;
% de=de(abs(de)<3);
title([num2str(mean(de),2) '±' num2str(std(de),ff) ', fit: ' num2str(fitp.b1,2) '±' num2str(ss,ff) ', in Gauss ' num2str(ingauss*100,'%2.0f') '%'])
xlim([-5*ss 5*ss])
end