classdef NPCLabelingEfficiency<interfaces.DialogProcessor&interfaces.SEProcessor
    properties
    end
    methods
        function obj=NPCLabelingEfficiency(varargin)        
            obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={};
            obj.showresults=true;
        end
        
        function out=run(obj,p)  
            out=runintern(obj,p);
            
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end

    end
end


function out=runintern(obj,p)
global rank lte rout
% p.bootstrap=1;
se=obj.SE;
quantifyevaluator='NPCLabelingQuantify_s';
% quantifyevaluator='NPCLabelingQuantify';
fields={'evaluation',quantifyevaluator};
fields2={'evaluation','generalStatistics'};

numbercornerassinged=getFieldAsVector(se.sites,fields{:},'numbercornerassigned');
numbercornerassigneddirect=getFieldAsVector(se.sites,fields{:},'numbercornerassigneddirect');
filenumber=getFieldAsVector(se.sites,'info','filenumber');
psf=getFieldAsVector(se.sites,fields2{:},'PSFlayers');

% fields3={'evaluation',quantifyevaluator,'timing'};
% timepoints=getFieldAsVectorInd(se.sites,fields3{:},'timepoints');
% nstart=getFieldAsVectorInd(se.sites,fields3{:},'nstart');
% nend=getFieldAsVectorInd(se.sites,fields3{:},'nend');
% frames=getFieldAsVectorInd(se.sites,fields{:},'coordinates','frame');

% nchunks=getFieldAsVectorInd(se.sites,fields3{:},'nchunks');
% dind=getFieldAsVectorInd(se.sites,fields3{:},'dind');

radius=getFieldAsVector(se.sites,fields{:},'radius');
siteid=getFieldAsVector(se.sites,'ID');
numlocs=getFieldAsVector(se.sites,fields{:},'numlocsf');
numlocsR=getFieldAsVector(se.sites,fields{:},'numlocsR');

use=getFieldAsVector(se.sites,'annotation','use');

if p.psfcheck
    indgood=psf>=p.PSFrange(1)&psf<=p.PSFrange(2);
else
    indgood=true(size(psf));
end

if p.filecheck
    indf=false(size(filenumber));
    for k=1:length(p.filenumbers)
        indf=indf | filenumber==p.filenumbers(k);
    end
    indgood=indgood&indf;
    filefile=p.filenumbers(1);
else
    filefile=se.sites(find(use,1,'first')).info.filenumber;
end
indgood=indgood&use;
indgood=indgood&~isnan(numbercornerassinged);

nb=0:p.corners;
numbercornerassinged=numbercornerassinged(indgood);
numbercornerassigneddirect=numbercornerassigneddirect(indgood);
radius=radius(indgood);
siteid=siteid(indgood);
numlocs=numlocs(indgood);
numlocsR=numlocsR(indgood);

% timepoints=timepoints(indgood,:);
% nstart=nstart(indgood,:);
% nend=nend(indgood,:);
% frames=frames(indgood,:);
psf=psf(indgood);

% nchunks=nchunks(indgood,:);
% dind=dind(indgood,:);
if isempty(numbercornerassinged) || length(numbercornerassinged)<=5
    warndlg('Not sufficient (>5) number of evaluations found. Make sure right evaluator used. Try redraw all.')
    return
end
if isfield(se.sites(1).evaluation.(quantifyevaluator),'numcornersfiltered_gt')
    gtexist=true;
    numcorners=getFieldAsVector(se.sites,fields{:},'numcornersfiltered_gt');
    numcornersa=getFieldAsVector(se.sites,fields{:},'numcornersall_gt');
%     fields3gt={'evaluation',quantifyevaluator,'timing_gt'};
%     timepoints_gt=getFieldAsVectorInd(se.sites,fields3gt{:},'timepoints');
%     nstart_gt=getFieldAsVectorInd(se.sites,fields3gt{:},'nstart');
%     nend_gt=getFieldAsVectorInd(se.sites,fields3gt{:},'nend');
%     timepoints_gt=timepoints_gt(indgood,:);
%     nstart_gt=nstart_gt(indgood,:);
%     nend_gt=nend_gt(indgood,:);
else
    gtexist=false;
end

ax0=obj.initaxis('Summary');
if any(psf>0)
axpsf=obj.initaxis('PSF');
histogram(axpsf,psf);
title(axpsf,['PSF range: ' num2str(p.PSFrange)])
xlabel('average PSF (nm)')
end

 ax3=obj.initaxis('corners');
ha=hist(numbercornerassigneddirect,nb);
p.ploton=true;
ph=p;
ph.ploton=false;
bs_assigned=bootstrp(20,@fitNPClabeling,numbercornerassigneddirect,ph);
berr_assigned=std(bs_assigned);
    bar(nb,ha)
    hold on
    p.ploton=true;
    pf=fitNPClabeling(ha,p);
    title(['assigned: ' num2str(100*pf,'%2.1f') '\pm' num2str(100*berr_assigned,'%2.1f')])
   axis tight
   results(1).assigned=pf*100; 
   results(2).assigned=berr_assigned*100;

    
% analyze time dependence
%     axtt=obj.initaxis('time');
%     [le,lte,nerre,rank]=evaluatetime(nend,p,'r');
%     hold(axtt,'on')
%     
%     clear nchall
%     tp=size(timepoints,2);
%     for k=1:tp-1
%         nh=nchunks(dind==k);
%         nchall{k+1}=nh(:);
%     end
%     [lc,ltc,nerrc,rank]=evaluatetime(nchall,p,'b');    
% 
%     title(axtt,['LE dye: ' num2str(le(3)*100,'%4.1f, ') ' , PAFP: ' num2str(lc(1)*100,'%4.1f, ')]);
%     legend('dye','model fit ' ,...
%         'pafp','lin fit ' ,'Location','southeast')
%     ylabel('labeling efficiency')
 
%     axvar=obj.initaxis('variance');
    %look at the variance when analyzing subset of sites only
%     nsmall=50; % size of subset
%     %first dye: nend + exp fit
%     ind=1:nsmall:size(nend,1)+1;
%     ph=p;
%     ph.ploton=false;
%     ph.bootstrap=false;
%     for k=length(ind)-1:-1:1
%         ledyeh=evaluatetime(nend(ind(k):ind(k+1)-1,:),ph,'g');
%         ledye(k)=ledyeh(3);
%     end
%     if isempty(rout)
%         indp=0;
%     else
%         indp=length(rout.dyemean);
%     end
%     rout.dyemean(indp+1)=mean(ledye);
%     rout.dyestd(indp+1)=std(ledye);
    
%     now chunks for pafp
    
%     tp=size(timepoints,2);
%     for k=length(ind)-1:-1:1
%         for l=1:tp-1
%             nchh=nchunks(ind(k):ind(k+1)-1,:);
%             dindh=dind(ind(k):ind(k+1)-1,:);
%             nh=nchh(dindh==l);
% 
%             nchallss{l+1}=nh(:);
%         end
%         
%         lepafph=evaluatetime(nchallss,ph,'k');
%         lepafp(k)=lepafph(1);
%         leall(k)=fitNPClabeling(nh,ph);
%     end
%     rout.pafpmean(indp+1)=mean(lepafp);
%     rout.pafpstd(indp+1)=std(lepafp); 
%     rout.bootstrapmean(indp+1)=pf;
%     rout.bootstraperr(indp+1)=berr_assigned;
%     rout.allmean(indp+1)=mean(leall);
%     rout.allstd(indp+1)=std(leall); 
    
%     hold(axvar,'off')
%     errorbar(axvar,rout.dyemean,rout.dyestd,'o')
%     hold(axvar,'on')
%     errorbar(axvar,rout.pafpmean,rout.pafpstd,'x')
    
    %number of localizations
    %aware of bias for not segmenting npcs with low number of localizations
axle=obj.initaxis('localizations');
hold off
hlocs=histogram(axle,numlocs,0:1:max(numlocs));
ns=floor(mean(numlocs)-1.2*std(numlocs));
% ns=1
fp=fit(hlocs.BinEdges(ns:end-1)',hlocs.Values(ns:end)','gauss1');
hold on
plot(hlocs.BinEdges+hlocs.BinWidth/2,fp(hlocs.BinEdges),'r')
plot(hlocs.BinEdges(ns:end-1)+hlocs.BinWidth/2,fp(hlocs.BinEdges(ns:end-1)),'r*')
title(axle,['mean: ' num2str(mean(numlocs),'%2.1f') ', LE:' num2str(mean(numlocs)/32*100,'%2.1f') ,' fit:' num2str(fp.b1,'%2.1f') ', LE:' num2str(fp.b1/32*100,'%2.1f')]);

%radii
axsize=obj.initaxis('size');
plotSElink(radius,numlocs,siteid,se,'o')
xlabel('radius (nm)')
ylabel('number of localizations')
results(3).assigned=0;
if gtexist %not from simulation
ax6=obj.initaxis('GT filtered');
%     p.ploton=false;
    bs_gtf=bootstrp(20,@fitNPClabeling,numcorners(indgood),ph);
    berr_gtf=std(bs_gtf);
    p.ploton=true;    
    
    hold off
    hnc=hist(numcorners(indgood),nb);
    bar(nb,hnc)
    hold on
    pf=fitNPClabeling(hnc,p);
    results(1).gtfilt=pf*100;
    results(2).gtfilt=berr_gtf*100;
    title(['GT filt: ' num2str(100*pf,'%2.1f') '\pm' num2str(100*berr_gtf,'%2.1f')])
    axis tight

%     ax6b=obj.initaxis('GT all');
%         p.ploton=false;
    bs_gtfa=bootstrp(20,@fitNPClabeling,numcornersa(indgood),ph);
    berr_gtfa=std(bs_gtfa);
    p.ploton=false;  
    
%     hold off
%     hnc=hist(numcornersa(indgood),nb);
%     bar(nb,hnc)
%     hold on
    pf=fitNPClabeling(hnc,p);
    results(1).gtall=pf*100; 
    results(2).gtall=berr_gtfa*100;
%     title(['GT all: ' num2str(100*pf,'%2.1f') '\pm' num2str(100*berr_gtfa,'%2.1f')])
%     axis tight
    numcorners=numcorners(indgood);
    rd=(0:length(numcorners)-1)/length(numcorners)/2;
    
%     ax7=obj.initaxis('comparison');
%     plot(numcorners+rd,numbercornerassigneddirect -numcorners,'d')
%     numbercornerassinedm=mean(numbercornerassinged-numcorners);
%     numbercornerassineddm=mean(numbercornerassigneddirect-numcorners);
%     title(['asssigned d: ' num2str(numbercornerassineddm,2)])  
    sp=obj.locData.files.file(filefile).info.simulationParameters;
    % results(3).le_gt=sp.labeling_efficiency;
    results(3).blinks_gt=sp.blinks;
    results(3).assigned=sp.photons;
    results(3).lifetime_gt=sp.lifetime;
    results(3).gtfilt=sp.background;

    results(3).gtall=sp.labeling_efficiency;
    results(1).lifetime_gt=0;results(2).lifetime_gt=0;
    results(1).blinks_gt=0;results(2).blinks_gt=0;
end

if p.filecheck
    results(1).file=(p.filenumbers(1));
    results(2).file=(p.filenumbers(end));
    results(3).file=0;

end

axp=ax0.Parent;
delete(ax0);
ht=uitable('Parent',axp);
struct2uitable(ht, results,'flip','%2.1f')
copytoexcel(results,'flip');
ht.ColumnWidth={50,40,30,40};
ht.ColumnName={};
out=[];

if p.copy2page
    sm=2;
    sn=2;
    f=figure;
    f.Renderer='painters';
    ht2=ht.copy;
    ht2.Parent=f;
    axtemp= subplot(sm,sn,[1 ]);
    ht2.Units='normalized';
    ht2.Position=axtemp.Position;
    delete(axtemp)
    if exist('axpsf','var')
        axpsf2=axpsf.copy;
        axpsf2.Parent=f;
        subplot(sm,sn,3,axpsf2)
    end
%     axt=ax1.copy;
%     axt.Parent=f;
%     subplot(sm,sn,4,axt)
%     axt=ax2.copy;
%     axt.Parent=f;
%     subplot(sm,sn,8,axt)
        axt=ax3.copy;
    axt.Parent=f;
    subplot(sm,sn,2,axt)
%     
%     axttime=axtt.copy;
%     axttime.Parent=f;
%     subplot(sm,sn,[5 6 9 10],axttime)   
    
    axt=axle.copy;
    axt.Parent=f;
    subplot(sm,sn,4,axt)   

    if exist('ax6','var')
        axt=ax6.copy;
    axt.Parent=f;
    subplot(sm,sn,5,axt) 
    end
    if exist('ax6b','var')
        axt=ax6b.copy;
    axt.Parent=f;
    subplot(sm,sn,6,axt) 

    end
    fn=se.files(p.filenumbers(1)).name;
    t=textwrap({fn},50);
%     th=text(axttime,0.15,0.05,t,'FontSize',9,'Interpreter','none');
end
filen=se.files(filefile).name;
%filename   LE  LEerrbs     numberofnpcs    locspernpcmean  locspernpcfit
   clipboard('copy',[filen sprintf(['\t' num2str(pf) '\t' num2str(berr_assigned) '\t' num2str(sum(indgood)) '\t' num2str(mean(numlocs)) '\t' num2str(fp.b1) ])])
   display(sprintf('filename  \t LE \t LEerrbs \t  numberofnpcs \t locspernpcmean \t locspernpcfit'))
end


function [out,pf,nerr,tp]=evaluatetime(n,pin,col)
if nargin <3
    col='k';
end
p=pin;
p.ploton=false;
p.fitrange=[1 8];
nerr=zeros(1,size(n,2))/100;
pf=zeros(1,size(n,2));
% fitstarts=[1 1 1 1 1 1 2 2 3];
for k=1:size(n,2)
    if iscell(n)
        nh=n{k};
    else
        nh=n(:,k);
    end
    
    if sum(nh>0)>5 %minium for fitting)
%         p.fitrange(1)=fitstarts(floor(mean(nh))+1);
        pf(k)=fitNPClabeling(nh,p);
        if p.bootstrap
            bs_numfoundint=bootstrp(10,@fitNPClabeling,nh,p);
            nerr(k)=std(bs_numfoundint);
        end
        nm(k)=mean(nh);
    else
        pf(k)=0;
        nerr(k)=0;
        nm(k)=0;
    end
end
tp=linspace(0,1,size(n,2)); %rank

% fitrange=pf<0.95;
% fitrange=nm>2 | pf==0;

% nm

if p.bootstrap
    delta=min(nerr(nerr>0));
    w=1./(nerr+delta);
%     w=w(fitrange);
    w=w/max(w);
else
    w=[];
end

%  fr=fit(tp(fitrange)',pf(fitrange)','poly1','Weights',w');
 

if pf(end)-pf(1) < 0 
%     fitrange=tp<=0.5;
    fitrange=nm>2;
    fitrange(end)=false;
    if ~isempty(w)
        w=w(fitrange);
    end
    [le, errle]=lscov(1-tp(fitrange)',pf(fitrange)',w);
%     fr=fit(tp(fitrange)',pf(fitrange)','poly1','Weights',w');
%     le=fr(0);
    %try other fit rank = a (pf(1-b ln(pf))
    g2 = fittype( @(le,a,x) logmodel(le,a,x));
%     fitrange(end)=false;
    fx2=fit(pf(fitrange)',1-tp(fitrange)',g2,'StartPoint',[pf(1) 0]);
%     findzero=pf(1)-0.05:0.001:pf(1)+0.05;
%     r0=fx2(findzero); 
%     ind=find(r0<0,1,'first');
%     lelog0=findzero(ind); 
%     lelog0
    lelog=fx2.le;
%     lelog
%     g1=fittype(@(le,a,b,x) le*x.*(1+a*x.*exp(-b*x)));
%     fx1=fit(1-tp(fitrange)',pf(fitrange)',g1,'StartPoint',[pf(1) 0.1 10]);
    fpol=fit(1-tp(fitrange)',pf(fitrange)','poly2','Weights',w');
    lepol=fpol(1);
    
%     figure(88);th=1-tp(fitrange)';pfh=pf(fitrange)';
%     hold off
% %     plot(th,pfh,'o',th,fx1(th)')
%     
%     plot(th,pfh,'x',fx2(pfh),pfh)
%     hold on
%     plot(th,fpol(th));
%     legend('dat','exp','log','pol')
%     ci=confint(fr);
%     dr=ci(2,:)-ci(1,:);
%     errle=sqrt(sum(dr.^2));


else
    lelog=0;
    lepol=0;
    fitrange=nm>2;
    if ~isempty(w)
        w=w(fitrange);
    end
    [le, errle]=lscov(tp(fitrange)',pf(fitrange)',w);
%      le=fr(tp(end));
end



if pin.ploton
    
    if pf(end)-pf(1) < 0 
        if all(nerr==0)
            plot(1-tp(fitrange),pf(fitrange),[col 'o'])
        else
            errorbar(1-tp(fitrange),pf(fitrange),nerr(fitrange),[col '.'])
        end
        hold on
        tpa=min(1-tp(fitrange)):0.01:1;
        leplot=0:0.01:fx2.le+0.02;
%         plot(tpa,fpol(tpa),col);
        plot(fx2(leplot),leplot,[col ]);
%         xlim([fx2f(pf(1)) 1])
    else
        if all(nerr==0)
            plot(tp(fitrange),pf(fitrange),[col 'x'])
        else
            errorbar(tp(fitrange),pf(fitrange),nerr(fitrange),[col '.'])
        end
        hold on
        plot(tp,le*(tp),col);
    end
xlabel('fraction used');
end
out=[le;errle;lepol];
out=[le;errle;lelog];
end



function out=logmodel(le,a,x)
    out=(x/le.*(1-a*log(x/le)));
    out(x==0)=0;

end
     

function pard=guidef(obj)
pard.t1.object=struct('String','Corners:','Style','text');
pard.t1.position=[1,1];

pard.corners.object=struct('String','8','Style','edit');
pard.corners.position=[1,2];
pard.corners.Width=0.5;

pard.t2.object=struct('String','Proteins/Corner','Style','text');
pard.t2.position=[1,3];

pard.rings.object=struct('String','4','Style','edit');
pard.rings.position=[1,4];
pard.rings.Width=0.5;


pard.psfcheck.object=struct('String','PSF range','Style','checkbox');
pard.psfcheck.position=[2,1];

pard.PSFrange.object=struct('String','80 150','Style','edit');
pard.PSFrange.position=[2,2];
pard.PSFrange.Width=1;

pard.t3.object=struct('String','fit range histogram','Style','text');
pard.t3.position=[2,3];

pard.fitrange.object=struct('String','3 8','Style','edit');
pard.fitrange.position=[2,4];
pard.fitrange.Width=1;

pard.filecheck.object=struct('String','filenumbers','Style','checkbox');
pard.filecheck.position=[3,1];

pard.filenumbers.object=struct('String','1:100','Style','edit');
pard.filenumbers.position=[3,2];
pard.filenumbers.Width=1;

% pard.bootstrap.object=struct('String','Bootstrap error bars for time analysis','Style','checkbox');
% pard.bootstrap.position=[4,1];
% pard.bootstrap.Width=2;

pard.copy2page.object=struct('String','Copy to own page','Style','checkbox');
pard.copy2page.position=[6,1];
pard.copy2page.Width=2;

pard.plugininfo.type='ROI_Analyze';

end