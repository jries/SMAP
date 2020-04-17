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
jaccard = truepositives / (truepositives + falsepositives + falsenegatives);
f1score = 2 * recall * precision / (recall + precision);

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


xbin=min(2,ceil(quantile(locerr.xerr,0.1)*20)/100);
ybin=min(2,ceil(quantile(locerr.yerr,0.1)*20)/100);
zbin=min(5,ceil(quantile(locerr.zerr,0.1)*20)/100);

subplot(3,4,1,'Parent',f)
[~, fitx] = fithistr(dx,xbin,[' ' num2str(mean(locerr.xerr),3)]);
xlabel('dx')
subplot(3,4,2,'Parent',f)
[~, fity] = fithistr(dy,ybin,[' ' num2str(mean(locerr.yerr),3)]);
xlabel('dy')
if isz
    subplot(3,4,3,'Parent',f)
    fitz=fithistr(dz,zbin,[' ' num2str(mean(locerr.zerr),3)]);
    xlabel('dz')
end
subplot(3,4,4,'Parent',f)
[~, fitphot] = fithistr(dphot,.005);
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
[fitxr, ~] = fithistr(dxr,0.1);
xlabel('dx/sqrt(CRLBx)')
subplot(3,4,6,'Parent',f)
[fityr, ~] = fithistr(dyr,0.1);
xlabel('dy/sqrt(CRLBy)')


subplot(3,4,8,'Parent',f)
fithistr(dphotr,.25);
xlabel('phot/sqrt(CRLBphot)')

if isz
    dzshift=(dz-fitz.b1);
    RMSEax=sqrt(mean(dzshift.^2));
    RMSEvol=sqrt(mean(dxshift.^2+dyshift.^2+dzshift.^2));
    dzr=(dz-fitz.b1)./locerr.zerr(inderr(indinz));
    subplot(3,4,7,'Parent',f)
[fitzr, ~] = fithistr(dzr,0.1);
xlabel('dz/sqrt(CRLBz)')



    subplot(3,4,10,'Parent',f)
    hold off

    q=quantile(vertcat(unique(locfit.bg(iBa)),unique(locgt.bg(iAa))),[0.02 0.98]);
    edges=floor(q(1)):1:ceil(q(2));
    [h1,edges1]=histcounts(locfit.bg(iBa),edges);
    [h2,edges2]=histcounts(locgt.bg(iAa),edges);
    h1(end+1)=0;
    h2(end+1)=0;
    
    bar(edges,h2/max(h2))

    hold on
    bb=bar(edges,h1/max(h1));
    bb.FaceAlpha=0.5;
else
    RMSEax=0;
    RMSEvol=0;
end
xlabel('bg')
% hold off
% dscatter(locgt.bg(iAa),locfit.bg(iBa))
% hold on
% plot([min(locgt.bg(iAa)),max(locgt.bg(iAa))],[min(locgt.bg(iAa)),max(locgt.bg(iAa))],'m')
% % ,'.',locgt.z(iAa),locgt.z(iAa),'k')
% xlabel('bg gt')
% ylabel('bg fit')

if ~all(locfit.z==0)
subplot(3,4,11,'Parent',f)
hold off
dscatter(locgt.z(iAa),locfit.z(iBa))
hold on
plot([min(locgt.z(iAa)),max(locgt.z(iAa))],[min(locgt.z(iAa)),max(locgt.z(iAa))],'m')
% ,'.',locgt.z(iAa),locgt.z(iAa),'k')
xlabel('z gt')
ylabel('z fit')

title(['RMSElat = ' num2str(RMSElat,3) ' nm, RMSEvol = ' num2str(RMSEvol,3) ' nm'])
end

subplot(3,4,12,'Parent',f);
if length(unique(locgt.phot(iAa)))>1
hold off
dscatter(locgt.phot(iAa),locfit.phot(iBa))
hold on
plot([min(locgt.phot(iAa)),max(locgt.phot(iAa))],[min(locgt.phot(iAa)),max(locgt.phot(iAa))],'m')
xlabel('phot gt')
ylabel('phpt fit')
ylim([0 1.7*max(locgt.phot(iAa))])
else
    histogram(locfit.phot(iBa))
end

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

% some higher level metrics
effcy_lat = efficiency(jaccard, RMSElat, 1.0);
effcy_ax = efficiency(jaccard, RMSEax, 0.5);
effcy_vol = (effcy_lat + effcy_ax) / 2;

sigma_crx = fitxr.c1 / sqrt(2);
sigma_cry = fityr.c1 / sqrt(2);
if isz
    sigma_crz = fitzr.c1 / sqrt(2);
    sigma_crvol = sqrt(sigma_crx ^ 2 + sigma_cry ^ 2 + sigma_crz ^ 2) / sqrt(3);  % is this correct?
else
    sigma_crz=0;
    sigma_crvol=0;
end

sigma_crlat = sqrt(sigma_crx ^ 2 + sigma_cry ^ 2) / sqrt(2);
sigma_crax = sigma_crz;


%add all results
results.precision = precision;
results.recall = recall;
results.jaccard = jaccard;
results.f1score = f1score;

results.truepositives = truepositives;
results.falsepositives = falsepositives;
results.falsenegatives = falsenegatives;

results.rmse_lat = RMSElat;
results.rmse_ax = RMSEax;
results.rmse_vol = RMSEvol;

results.effcy_lat = effcy_lat;
results.effcy_ax = effcy_ax;
results.effcy_vol = effcy_vol;

% add the err / crlb stuff here
results.dx_cr = sigma_crx;
results.dy_cr = sigma_cry;
results.dz_cr = sigma_crz;

results.lat_cr = sigma_crlat;
results.ax_cr = sigma_crax;
results.vol_cr = sigma_crvol;

end

function e = efficiency(jac, rmse, alpha)
    %%%
    % Calculate efficiency following Sage et al. 2019, superres fight club.
    % Alpha_lat = 1, Alpha_ax = 0.5nm
    % jac should be in 0...1 not percent.
    % :return efficiency in 0...1 range.
    %%%
    e = (100 - ((100 * (1 - jac)) ^ 2 + alpha ^ 2 * rmse ^ 2) ^ 0.5) / 100;
    return
end

function [fitp, fitp2] = fithistr(de,dn,addtxt)
if nargin<3
    addtxt='';
end
ff='%1.1f';
ff2='%1.2f';
qq=quantile(de,[0.002 0.998]);
if qq(1)==qq(2) || any(isinf(abs(qq)))
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
n2=(max(1,mp-dn2):min(mp+dn2,length(hn)))';

% fitp2=fit(nf(n2)',hn(n2)','gauss2','StartPoint',[fitp.a1,fitp.b1,fitp.c1,fitp.a1/10,fitp.b1,fitp.c1*10]);
if numel(n2) > 3
    fitp2=fit(double(nf(n2)'),double(hn(n2)'),'gauss1','StartPoint',[fitp.a1,fitp.b1,fitp.c1]);
else
    fitp2 = [];
end

hold on
plot(nf,fitp(nf),'g')
if ~isempty(fitp2)
plot(nf(n2),fitp2(nf(n2)),'r')
end

fituse=fitp;
ss2=fituse.c1/sqrt(2);
ingauss=fituse.a1*sqrt(pi)*fituse.c1/length(de)/dn;
% de=de(abs(de)<3);
title([num2str(mean(de),2) '±' num2str(std(de),ff) ', fit: ' num2str(fituse.b1,ff) '±' num2str(ss2,ff2) ', in Gauss ' num2str(ingauss*100,'%2.0f') '%' addtxt])
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