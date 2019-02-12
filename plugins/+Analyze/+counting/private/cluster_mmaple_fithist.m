function parout=cluster_mmaple_fithist(par,histogram)

if nargin==0
    par.N0_v=10;
    par.pmature_v=0.5;
    par.pblink_v=0.2;
    par.monomer_v=0.0;
    par.ncluster_v=0;
    par.blinkmode.selection='Exponential';
    par.lenbthscale=15;
    par.N0_fit=0;
    par.pmature_fit=1;
    par.pblink_fit=1;
    par.monomer_fit=0;
    par.ncluster_fit=0;
    par.bin=1;
    
end
% %evaluate cluster mmaple

axhist=initaxis(par.resultstabgroup,'histogram');
axwhist=initaxis(par.resultstabgroup,'weighted histogram');
axcum=initaxis(par.resultstabgroup,'cumulative');

x=histogram.c;
maxfit=min(par.fitrange_max,length(x));
minfit=max(1,par.fitrange_min);

% hnorm=histogram.h/sum(histogram.h(:));

histogram.cum=cumsum(histogram.h(minfit:maxfit));
% histogram.cum=histogram.cum/histogram.cum(end);
histogram.cumx=histogram.c(minfit:maxfit);
switch par.fitselection.Value
    case 1 %cum
        y=histogram.cum;
        fitfun=@cumforfit;
        fitamplitude=1;
        fitax=axcum;
        x=histogram.cumx;
        maxfit=maxfit-minfit+1;
        minfit=1;
        
     
    case 2 %hist
        y=histogram.h;
        fitfun=@histforfit;
        fitamplitude=1;
        fitax=axhist;
    case 3 %whist
        y=sqrt(histogram.h);
        fitfun=@whistforfit;
         fitamplitude=1;
         fitax=axwhist;
end

ampstart=sum(histogram.h(:));
% x=histogram.c;
% y=histogram.h;
% y=cumsum(histogram.h);
% y=y/y(end);

% maxfit=min(par.fitrange_max,length(x));
% minfit=max(1,par.fitrange_min);

% dx=x(2)-x(1);
%fit theoretical distribution
% ft = fittype('a*brightnessdistribution( x, 10,1,0 )');
% fit(x,y,ft)

allpar=[par.N0_v,par.pmature_v,par.pblink_v,par.ncluster_v,ampstart];
fitind=[par.N0_fit,par.pmature_fit,par.pblink_fit,par.ncluster_fit,fitamplitude];

%amplitude not fitted: converges faster. But does not exlude outliers. Maybe adjust later?

% fitp=lsqcurvefit(@bforfit,startp,x(minfit:maxfit)',y(minfit:maxfit)',[0 0 ],[inf,inf]);
addpar.blinkmode=par.blinkmode.selection;



if par.blinkmode.Value==1
    ubblink=inf;
else
    ubblink=1;
end
fitp=fitSome(fitfun,x(minfit:maxfit)',y(minfit:maxfit)',allpar,fitind,[0 0.01 0.001 0 0],[inf,1, ubblink,2.5, inf],addpar);
% fitp=lsqcurvefit(@bforfit,startp,x(minfit:maxfit)',y(minfit:maxfit)',[0 0.01 0.001 0],[inf,1, 1,2.5]);


%plot histogram
dx=histogram.c(2)-histogram.c(1);
binning=par.bin;
hold(axhist,'off')
plothist(axhist,histogram.c,histogram.h, binning, @bar,'b');
% plothist(axhist,histogram.c,histogram.h, binning, @stairs,'b');
hold(axhist,'on')
plothist(axhist,histogram.c,histforfit(fitp,histogram.c',addpar)*dx, binning,[], 'r');
plothist(axhist,histogram.c,histforfit(allpar,histogram.c',addpar)*dx, binning, [], 'k--');

hold(axwhist,'off')
plothist(axwhist,histogram.c,sqrt(histogram.h), binning, @stairs,'b');
hold(axwhist,'on')
plothist(axwhist,histogram.c,whistforfit(fitp,histogram.c',addpar), binning, [],'r');
plothist(axwhist,histogram.c,whistforfit(allpar,histogram.c',addpar), binning, [], 'k--');

hold(axcum,'off')
plothist(axcum,histogram.cumx,histogram.cum, binning,@stairs, 'b');
hold(axcum,'on')
plothist(axcum,histogram.cumx,cumforfit(fitp,histogram.cumx',addpar), binning, [],'r');
plothist(axcum,histogram.cumx,cumforfit(allpar,histogram.cumx',addpar), binning,[],  'k--');

% 
% plot(histogram.c,histogram.h,'b',histogram.c,histforfit(fitp,histogram.c',addpar),'r',histogram.c,histforfit(allpar,histogram.c',addpar),'k--','Parent',axhist);
% plot(histogram.c,sqrt(histogram.h),'b',histogram.c,whistforfit(fitp,histogram.c',addpar),'r',histogram.c,whistforfit(allpar,histogram.c',addpar),'k--','Parent',axwhist);
% plot(histogram.cumx,histogram.cum,'b',histogram.cumx,cumforfit(fitp,histogram.cumx',addpar),'r',histogram.cumx,cumforfit(allpar,histogram.cumx',addpar),'k--','Parent',axcum);
% 
% fitax.NextPlot='add';
% plot(x,fitfun(allpar,x',addpar),'k--','Parent',fitax);
% fitax.NextPlot='replace';
tit{1}=num2str(fitp,3);
meanlocs=sum(histogram.h.*histogram.c)/sum(histogram.h);
hf=histforfit(fitp,histogram.c',addpar);
meanlocsfit=sum(hf.*histogram.c')/sum(hf);
tit{2}=['mean locs: ' num2str(meanlocs,'%3.1f') ' (fit: ' num2str(meanlocsfit,'%3.1f') '), locs/protein: ' num2str(meanlocs/fitp(1),'%2.2f') ' (fit: ' num2str(meanlocsfit/fitp(1),'%2.2f') ')'];
axhist.Title.String=tit;
axcum.Title.String=tit;
axwhist.Title.String=tit;
initaxis(fitax.Parent,[],'keep');

% 
% plot(x,fitfun(fitp,x',addpar),'r','Parent',fitax)
% 
% hold off
% title(fitp);
% 
% plot(x,y,'Parent',fitax)
% hold on
% plot(x,fitfun(allpar,x',addpar),'b--','Parent',fitax)
% plot(x,histogram.h)
% hold on
% plot(x,histforfit(allpar,x',addpar),'b--')
% plot(x,histforfit(fitp,x',addpar),'r')
% hold off
%assign numlocsreduced to each localization

% hall(hindex).histogram=histogram;
% hall(hindex).histfit.y=bforfit(fitp,x',addpar);
% hall(hindex).fitp=fitp;
% hall(hindex).cluster=cluster;

% plotc.x=vertcat(cluster(:).x);
% plotc.y=vertcat(cluster(:).y);
% plotc.c=vertcat(cluster(:).locsdfv);
% plotc.s=vertcat(cluster(:).locprec);

parout=par;
parout.pmature_v=fitp(2);
parout.pblink_v=fitp(3);
% parout.monomer_v=0.2;
parout.ncluster_v=fitp(4);
parout.meanlocs = num2str(meanlocs,'%3.1f'); % added by Yu-le 
% recgui.initaxis(par.resultstabgroup,'cumulative','keep');

if par.N0_fit
    initaxis(par.resultstabgroup,'N0fit');
N0range=max(1,par.N0_v-5):par.N0_v+5;
% options= optimoptions('lsqcurvefit');
% resnorm=0*N0range;

hold off
plot(x,y,x,fitfun(fitp,x',addpar))

res=zeros(length(N0range),1);
for k=1:length(N0range)
    allpar=fitp;
    allpar(1)=N0range(k);
    fitind=[0 0 0 0  1];
    fitp2=fitSome(fitfun,x(minfit:maxfit)',y(minfit:maxfit)',allpar,fitind,[0 0.01 0.001 0 0],[inf,1, ubblink,2.5, inf],addpar);
%     fitp(1)=N0range(k);
%     fitp(5)=1.7
    bN=fitfun(fitp2,x',addpar);
%     ratio=mean(y'./bN)
%     ratio=y'\bN;
%        bN=bN./ratio;
    hold on
    plot(x,bN)
    
    res(k)=sum((bN-y').^2);
%    startp=[1] ;
%    [fitp,resnorm(k)]=lsqcurvefit(@bforfitN0,startp,x(minfit:maxfit)',y(minfit:maxfit)',[0 ],[inf],options,N0range(k));
%    fita{k}=fitp;
%    disp([resnorm, N0range(k)])
end
initaxis(par.resultstabgroup,'N0fit:res');
plot(N0range,res);
[~,ind]=min(res);
parout.N0_v=N0range(ind);


% 
% [~,mink]=min(resnorm);
% figure(87)
% subplot(2,1,1)
% plot(N0range,resnorm);
% subplot(2,1,2)
% plot(x,y)
% hold on
% plot(x,bforfitN0(fita{mink},x',N0range(mink)),'r')
% hold off



% posa.c=posa.numlocsred;
% if 0
%     pixrec2=10;
% lut=jet;
% lut(1,:)=[0 0 0];
% pixrec2=5;
% % imsr=gaussrender(plotc,[0 cluster(1).xmax],[0 cluster(1).ymax],pixrec2, pixrec2,lut,[0 40]);
% figure(32)
% image(imsr)
% colormap jet
% colorbar
end
% axes(axcum);

function plothist(axhist,x,h, binning, plotfun,varargin)
if isempty(plotfun)
    plotfun=@plot;
end

    xb=x(1:binning:end);
    hb=bindata(x,h,xb);
     
    plotfun(xb,hb,varargin{:},'Parent',axhist);

% profile viewer

function fitpout=fitSome(fun,x,y,allpar,fitind,lbi,ubi,addpar)
if nargin<6
    lb=[];
    ub=[];
else
    lb=lbi(fitind==1);
    ub=ubi(fitind==1);
end
startp=allpar(fitind==1);

addp.fitind=fitind==1;
addp.allpar=allpar;
addp.fitfunction=fun;
addp.addpar=addpar;
options= optimoptions('lsqcurvefit');
options.Display='off';
fitp=lsqcurvefit(@callfitfun,startp,x,y,lb,ub,options,addp);
fitpout=allpar;
fitpout(fitind==1)=fitp;

function y=callfitfun(fitpar,x,addp)
par=addp.allpar;
par(addp.fitind)=fitpar;
addpar=addp.addpar;
y=addp.fitfunction(par,x,addpar);

% function y=bforfit(par,x,addpar)
% % par(5)=1;
% xi=0:x(end);
% yi=par(5)*brightnessdistribution(xi',round(par(1)),par(2),par(3),par((4)),'blinkmode',addpar.blinkmode);
% yii=cumsum(yi);
% y=yii(x(1)+1:end);
% 
% function y=bforfith(par,x,addpar)
% % par(5)=1;
% % xi=0:x(end);
% y=par(5)*brightnessdistribution(x,round(par(1)),par(2),par(3),par((4)),'blinkmode',addpar.blinkmode);
% % yii=cumsum(yi);
% % y=yii(x(1)+1:end);

function y=histforfit(par,x,addpar)
y=par(5)*brightnessdistribution(x,round(par(1)),par(2),par(3),par((4)),'blinkmode',addpar.blinkmode);

function y=whistforfit(par,x,addpar)
y=sqrt(par(5)*brightnessdistribution(x,round(par(1)),par(2),par(3),par((4)),'blinkmode',addpar.blinkmode));

function y=cumforfit(par,x,addpar)
% xi=0:x(end);
% xi=x';
% yi=par(5)*brightnessdistribution(x,round(par(1)),par(2),par(3),par((4)),'blinkmode',addpar.blinkmode);
yi=histforfit(par,x,addpar);
% yi=yi/sum(yi);
yii=cumsum(yi);
% y=yii(x(1)+1:end);
y=yii;
% function y=bforfitN0(1par,x,addpar)
% y=par(5)*brightnessdistribution(x,round(par(1)),par(2),par(3),par((4)),'blinkmode',addpar.blinkmode);
% y=cumsum(y);
% y=par(1)*brightnessdistribution(x,round(par(2)),.50,.33,0,'blinkmode','Exponential');
