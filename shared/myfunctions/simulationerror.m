function results=simulationerror(locgt,locfit,whicherr,searchradius)

f=gcf;
isz=true;
% pixelsize=100
if ~isfield(locgt,'x'), locgt.x=locgt.xnm; end
if ~isfield(locfit,'x'), locfit.x=locfit.xnm; end
if ~isfield(locgt,'y'), locgt.y=locgt.ynm; end
if ~isfield(locfit,'y'), locfit.y=locfit.ynm; end
if ~isfield(locgt,'z')
    if isfield(locgt,'znm')
        locgt.z=locgt.znm; 
    else
        locgt.z=0*locgt.x;
        isz=false;
    end
end
if ~isfield(locfit,'z')
    if isfield(locfit,'znm')
        locfit.z=locfit.znm; 
    else
        locfit.z=0*locfit.x;
        isz=false;
    end
end

if ~isfield(locgt,'N'), locgt.N=locgt.phot; end
if ~isfield(locfit,'N'), locfit.N=locfit.phot; end



if nargin<4 || isempty(searchradius)
    searchradius=200;
end

[iAa,iBa,nA,nB,nseen]=matchlocsall(locgt,locfit,0,0,searchradius(1));
totallocs=length(locgt.x);
falsepositives=length(nB);
falsenegatives=length(nA);
truepositives=length(iAa);
precision=truepositives/(truepositives+falsepositives);
recall=truepositives/(truepositives+falsenegatives);

matched=length(iAa);
if isz
dz=locgt.z(iAa)-locfit.z(iBa);
indinz=abs(dz)<searchradius(end);
dz=dz(indinz);
else
    indinz=true(size(iAa));
end
dx=locgt.x(iAa(indinz))-locfit.x(iBa(indinz));
dy=locgt.y(iAa(indinz))-locfit.y(iBa(indinz));

dphot=locgt.N(iAa(indinz))./locfit.N(iBa(indinz));

% if nargin<3 || isempty(psf)
%     % test if errors are transferred
%     [lp,errphot]=Mortensen(locgt.phot(iAa(indinz)),locgt.bg(iAa(indinz)),150, pixelsize,0);
% %     crlb=1000./sqrt(locgt.phot(iAa(indinz))); %
%     dxr=dx./lp; %xeems to be closer to 1
%     dyr=dy./lp;
%     dzr=dz./lp/3; 
% else    
%     crlb=psf.crlb(locgt.phot(iAa(indinz)),locgt.bg(iAa(indinz)),-locgt.z(iAa(indinz)));
%     dxr=dx./sqrt(crlb(:,2))/pixelsize; %xeems to be closer to 1
%     dyr=dy./sqrt(crlb(:,1))/pixelsize;
%     dzr=dz./sqrt(crlb(:,5));
% end

%errors for normalization
switch whicherr 
    case 1
        locerr=locgt;
        inderr=iAa;
    case 2
        locerr=locfit;
        inderr=iBa;
    otherwise
        disp('third argument should be 1 if error is taken from first argument, 2 if error is in second')
end
if ~isfield(locerr,'xerr')
    if isfield(locerr,'xnmerr')&&~isempty(locerr.xnmerr)
        locerr.xerr=locerr.xnmerr;
    elseif isfield(locerr,'locprecnm')&&~isempty(locerr.locprecnm)
        locerr.xerr=locerr.locprecnm;
    else
        locerr.xerr=Mortensen(locgt.phot,locgt.bg,150,100,0);
        disp('error in x estimated using Mortensen');
    end
end
if ~isfield(locerr,'yerr')
    if isfield(locerr,'ynmerr')&&~isempty(locerr.ynmerr)
        locerr.yerr=locerr.ynmerr;
    elseif isfield(locerr,'locprecnm')&&~isempty(locerr.locprecnm)
        locerr.yerr=locerr.locprecnm;
    else
        locerr.yerr=Mortensen(locgt.phot,locgt.bg,150,100,0);
        disp('error in y estimated using Mortensen');
    end
end
if ~isfield(locerr,'zerr') && isz
    if isfield(locerr,'znmerr')&&~isempty(locerr.znmerr)
        locerr.zerr=locerr.znmerr;
    elseif isfield(locerr,'locprecznm')&&~isempty(locerr.locprecznm)
        locerr.zerr=locerr.locprecznm;
    else
        locerr.zerr=Mortensen(locgt.phot,locgt.bg,150,100,0)*3;
        disp('error in z estimated using Mortensen');
    end
end


if ~isfield(locerr,'Nerr')
    if isfield(locerr,'photerr')&&~isempty(locerr.photerr)
        locerr.Nerr=locerr.photerr;
    else
        [~,locerr.Nerr]=Mortensen(locgt.phot,locgt.bg,150,100,0);
        
        disp('error in phtons estimated using Rieger');
    end
end


subplot(3,4,1,'Parent',f)
fitx=fithistr(dx,2);
xlabel('dx')
subplot(3,4,2,'Parent',f)
fity=fithistr(dy,2);
xlabel('dy')
if isz
    subplot(3,4,3,'Parent',f)
    fitz=fithistr(dz,5);
    xlabel('dz')
end
subplot(3,4,4,'Parent',f)
fitphot=fithistr(dphot,.01);
xlabel('dphot')
% shifted dx
dxshift=(dx-fitx.b1);
dyshift=(dy-fity.b1);

RMSElat=sqrt(mean(dxshift.^2+dyshift.^2));
%renormalized
dxr=(dx-fitx.b1)./locerr.xerr(inderr(indinz));
dyr=(dy-fity.b1)./locerr.yerr(inderr(indinz));

dphotr=(locgt.N(iAa(indinz))/fitphot.b1-locfit.N(iBa(indinz)))./locerr.Nerr(inderr(indinz));
subplot(3,4,5,'Parent',f)
fithistr(dxr,0.25);
xlabel('dx/sqrt(CRLBx)')
subplot(3,4,6,'Parent',f)
fithistr(dyr,0.25);
xlabel('dy/sqrt(CRLBy)')


subplot(3,4,8,'Parent',f)
fithistr(dphotr,.5);
xlabel('phot/sqrt(CRLBphot)')

if isz
    dzshift=(dz-fitz.b1);
    RMSEvol=sqrt(mean(dxshift.^2+dyshift.^2+dzshift.^2));
    dzr=(dz-fitz.b1)./locerr.zerr(inderr(indinz));
    subplot(3,4,7,'Parent',f)
fithistr(dzr,0.25);
xlabel('dz/sqrt(CRLBz)')
end

q=quantile(vertcat(unique(locfit.bg(iBa)),unique(locgt.bg(iAa))),[0.02 0.98]);
edges=floor(q(1)):1:ceil(q(2));
[h1,edges1]=histcounts(locfit.bg(iBa),edges);
[h2,edges2]=histcounts(locgt.bg(iAa),edges);
h1(end+1)=0;
h2(end+1)=0;
subplot(3,4,10,'Parent',f)
hold off
bar(edges,h2/max(h2))

hold on
bb=bar(edges,h1/max(h1));
bb.FaceAlpha=0.5;
xlabel('bg')
% hold off
% dscatter(locgt.bg(iAa),locfit.bg(iBa))
% hold on
% plot([min(locgt.bg(iAa)),max(locgt.bg(iAa))],[min(locgt.bg(iAa)),max(locgt.bg(iAa))],'m')
% % ,'.',locgt.z(iAa),locgt.z(iAa),'k')
% xlabel('bg gt')
% ylabel('bg fit')

subplot(3,4,11,'Parent',f)
hold off
dscatter(locgt.z(iAa),locfit.z(iBa))
hold on
plot([min(locgt.z(iAa)),max(locgt.z(iAa))],[min(locgt.z(iAa)),max(locgt.z(iAa))],'m')
% ,'.',locgt.z(iAa),locgt.z(iAa),'k')
xlabel('z gt')
ylabel('z fit')

title(['RMSElat = ' num2str(RMSElat,3) ' nm, RMSEvol = ' num2str(RMSEvol,3) ' nm'])

subplot(3,4,12,'Parent',f);

hold off
dscatter(locgt.phot(iAa),locfit.phot(iBa))
hold on
plot([min(locgt.phot(iAa)),max(locgt.phot(iAa))],[min(locgt.phot(iAa)),max(locgt.phot(iAa))],'m')
xlabel('phot gt')
ylabel('phpt fit')
ylim([0 1.7*max(locgt.phot(iAa))])

%matchlocs in defined region
% false pos, false neg
% matched: histogram dx/crlbx dy/crlby dz/crlbz
% std of these quantities vs phot, vs z
subplot(3,4,9,'Parent',f);
hold off
plot(locfit.x,locfit.y,'.',locgt.x,locgt.y,'.')
ff='%2.0f';
title(['FP: ' num2str(falsepositives/totallocs*100,ff), '%, FN: ' num2str(falsenegatives/totallocs*100,ff) '%'...
    ', P: ' num2str(precision,2) ', R: ' num2str(recall,2)]);
results.precision=precision;
results.recall=recall;
results.falsepositives=falsepositives;
%add all results

end

function fitp2=fithistr(de,dn)
ff='%1.1f';
ff2='%1.2f';
qq=quantile(de,[0.002 0.998]);
if qq(1)==qq(2)
    fitp2=zeros(5);
    return
end
n=floor(qq(1)):dn:ceil(qq(end));
hold off
histogram(de,n);
hn=histcounts(de,n);
nf=n(1:end-1)+(n(2)-n(1))/2;
fitp=fit(double(nf'),double(hn'),'gauss1');
ss=fitp.c1/sqrt(2);

% try fitting with second gauss
mp=find(nf>fitp.b1,1,'first');
dn2=min(mp-1,ceil(2*ss/dn));
n2=(mp-dn2:mp+dn2)';

% fitp2=fit(nf(n2)',hn(n2)','gauss2','StartPoint',[fitp.a1,fitp.b1,fitp.c1,fitp.a1/10,fitp.b1,fitp.c1*10]);
fitp2=fit(double(nf(n2)'),double(hn(n2)'),'gauss1','StartPoint',[fitp.a1,fitp.b1,fitp.c1]);

hold on
plot(nf,fitp(nf),'g')
plot(nf(n2),fitp2(nf(n2)),'r')

fituse=fitp;
ss2=fituse.c1/sqrt(2);
ingauss=fituse.a1*sqrt(pi)*fituse.c1/length(de)/dn;
% de=de(abs(de)<3);
title([num2str(mean(de),2) '±' num2str(std(de),ff) ', fit: ' num2str(fituse.b1,ff) '±' num2str(ss2,ff2) ', in Gauss ' num2str(ingauss*100,'%2.0f') '%'])
xlim([fituse.b1-5*ss fituse.b1+5*ss])
axh=gca;

wx=mean(axh.XLim(:)); dx=axh.XLim(2)-axh.XLim(1);
text(wx+dx/7,max(hn)*0.9,[ '\sigma=' num2str(ss2,ff2)],'FontSize',16)
t2=text(wx+dx/4,max(hn)*0.8,[num2str(ingauss*100,'%2.0f') '%'],'FontSize',16);
end


function [lp,errphot]=Mortensen(N,Bg,PSF, pixel,cmosn)
b=sqrt(Bg+cmosn^2);

PSFa=sqrt(PSF^2+pixel^2/12);
v=PSFa^2./N.*(16/9+8*pi*PSFa^2*b.^2./N/pixel^2);
lp=sqrt(v);


s_a=PSF/pixel; %sigmapsf/pixelsize
tau=2*pi*(Bg)*(s_a^2+1/12)./N;
errphot2=N.*(1+4*tau+sqrt(tau./(14*(1+2*tau)))); %This is Rieger...
errphot=sqrt(errphot2);
end