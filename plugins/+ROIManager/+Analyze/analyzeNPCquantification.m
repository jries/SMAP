classdef analyzeNPCquantification<interfaces.DialogProcessor&interfaces.SEProcessor
    properties
    end
    methods
        function obj=analyzeNPCquantification(varargin)        
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

% function pf=fithist(hi,n,p)
% 
% shi=sum(hi);
% % hi=hi/shi;
% 
% corners=p.corners;rings=p.rings;
% n1=find(n>=p.fitrange(1),1,'first');
% n2=find(n<=p.fitrange(2),1,'last');
% range=n1:n2;
% % x=(0:corners)';
% x=n(range)';
% % clusterfromlabeling(x,corners,rings,.5)
% ft=fittype('a*clusterfromlabeling(x,corners,rings,p)','problem',{'corners','rings'});
% f=fit(x,hi(range)',ft,'problem',{corners, rings},'Lower',[0 0.01],'Upper',[inf .99],'Start',[shi .4]);
% plot(n,f(n),'-g')
% plot(x,f(x),'-*r')
% pf=f.p;
% end

function out=runintern(obj,p)
% p.bootstrap=1;
se=obj.SE;
fields={'evaluation','NPCLabelingQuantify'};
fields2={'evaluation','generalStatistics'};

numfoundint=getFieldAsVector(se.sites,fields{:},'numfoundint');
numfoundrat=getFieldAsVector(se.sites,fields{:},'numfoundrat');
numbercornerassinged=getFieldAsVector(se.sites,fields{:},'numbercornerassigned');
numbercornerassigneddirect=getFieldAsVector(se.sites,fields{:},'numbercornerassigneddirect');
filenumber=getFieldAsVector(se.sites,'info','filenumber');
psf=getFieldAsVector(se.sites,fields2{:},'PSFlayers');
% psf=ones(1,length(psf)); %XXXXX

fields3={'evaluation','NPCLabelingQuantify','timing'};
timepoints=getFieldAsVectorInd(se.sites,fields3{:},'timepoints');
nstart=getFieldAsVectorInd(se.sites,fields3{:},'nstart');
nend=getFieldAsVectorInd(se.sites,fields3{:},'nend');
frames=getFieldAsVectorInd(se.sites,fields{:},'coordinates','frame');

radius=getFieldAsVector(se.sites,fields{:},'radius');
siteid=getFieldAsVector(se.sites,'ID');
numlocs=getFieldAsVector(se.sites,fields{:},'numlocsf');
numlocsR=getFieldAsVector(se.sites,fields{:},'numlocsR');
% numfoundint=numfoundint(use);
% numfoundrat=numfoundrat(use);
% numbercornerassinged=numbercornerassinged(use);
% psf=psf(use);
% 
% 
% 
% timepoints=timepoints(use,:);
% nstart=nstart(use,:);
% nend=nend(use,:);
% 
% timepoints_gt=timepoints_gt(use,:);
% nstart_gt=nstart_gt(use,:);
% nend_gt=nend_gt(use,:);

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

numfoundint=numfoundint(indgood);
numfoundrat=numfoundrat(indgood);
numbercornerassinged=numbercornerassinged(indgood);
numbercornerassigneddirect=numbercornerassigneddirect(indgood);
radius=radius(indgood);
siteid=siteid(indgood);
numlocs=numlocs(indgood);
numlocsR=numlocsR(indgood);

timepoints=timepoints(indgood,:);
nstart=nstart(indgood,:);
nend=nend(indgood,:);
frames=frames(indgood,:);

if isempty(numbercornerassinged) || length(numbercornerassinged)<=5
    warndlg('Not sufficient (>5) number of evaluations found. Make sure right evaluator used. Try redraw all.')
    return
end
if isfield(se.sites(1).evaluation.NPCLabelingQuantify,'numcornersfiltered_gt')
    gtexist=true;
    numcorners=getFieldAsVector(se.sites,fields{:},'numcornersfiltered_gt');
    numcornersa=getFieldAsVector(se.sites,fields{:},'numcornersall_gt');
    fields3gt={'evaluation','NPCLabelingQuantify','timing_gt'};
    timepoints_gt=getFieldAsVectorInd(se.sites,fields3gt{:},'timepoints');
    nstart_gt=getFieldAsVectorInd(se.sites,fields3gt{:},'nstart');
    nend_gt=getFieldAsVectorInd(se.sites,fields3gt{:},'nend');
    timepoints_gt=timepoints_gt(indgood,:);
    nstart_gt=nstart_gt(indgood,:);
    nend_gt=nend_gt(indgood,:);
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

ax1=obj.initaxis('from gap');
ax2=axes(ax1.Parent);
subplot(1,2,1,ax1);
    p.ploton=false;
    bs_numfoundint=bootstrp(20,@fitNPClabeling,numfoundint,p);
    berr_numfoundint=std(bs_numfoundint);
    p.ploton=true;    
hi=hist(numfoundint,nb);
    bar(nb,hi)
    hold on
    pf=fitNPClabeling(hi,p);
    title(['gap int: ' num2str(100*pf,'%2.1f') '\pm' num2str(100*berr_numfoundint,'%2.1f')])
    axis tight
    results.gapint=pf*100;   
    results(2).gapint=berr_numfoundint*100;
subplot(1,2,2,ax2);
    p.ploton=false;
    bs_numfoundrat=bootstrp(20,@fitNPClabeling,numfoundrat,p);
    berr_numfoundrat=std(bs_numfoundrat);
    p.ploton=true;    
    
    hold off
    hr=hist(numfoundrat,nb);
    bar(nb,hr)
    hold on
    pf=fitNPClabeling(hr,p);
    title(['gap frac: ' num2str(100*pf,'%2.1f')  '\pm' num2str(100*berr_numfoundrat,'%2.1f')])
    axis tight
    results(1).gapf=pf*100 ; 
    results(2).gapf=berr_numfoundrat*100;
    
 ax3=obj.initaxis('ass.+all');

ax4=axes(ax3.Parent);

   subplot(1,2,1,ax3);
    ha=hist(numbercornerassigneddirect,nb);
    p.ploton=false;
    bs_assigned=bootstrp(20,@fitNPClabeling,numbercornerassigneddirect,p);
    berr_assigned=std(bs_assigned);
%     [ncbc,bstat]=bootci(30,@fitNPClabeling,numbercornerassigneddirect,p);
%     ncbc=bootci(250,@fitNPClabeling,numbercornerassigneddirect,p);
    
    bar(nb,ha)
    hold on
    p.ploton=true;
    pf=fitNPClabeling(ha,p);
    title(['assigned: ' num2str(100*pf,'%2.1f') '\pm' num2str(100*berr_assigned,'%2.1f')])
   axis tight
   results(1).assigned=pf*100; 
   results(2).assigned=berr_assigned*100;
   
subplot(1,2,2,ax4);
    hall=hi+hr+ha;
    bar(nb,hall)
    hold on
    pf=fitNPClabeling(hall,p);
    title(['all: ' num2str(100*pf,'%2.1f')])
    axis tight 
    results(1).all=pf*100; 

    
% analyze time dependence
    %at some time we have to rescale frames. Already do so in Quantify? Or
    %only here? Time axis should be localizations, not frames, to get rid
    %of activation history
%     frames=frames(frames>0);
%     qq=timepoints(1,:)/timepoints(1,end);
%     timepointsq=myquantile(frames,qq);
    
    axtt=obj.initaxis('time');
    hold(axtt,'off')
    [ls,lts,nerrs,rank]=evaluatetime(nstart,p,'r');
    hold(axtt,'on')
    [le,lte,nerre,rank]=evaluatetime(nend,p,'b');
    
    nerra=sqrt(nerrs.^2+nerre.^2);
    nerra(nerra==0)=min(nerra(nerra>0));
    lta=lts+lte;
    mlta=sum(lta./nerra)/sum(1./nerra);
    
    plot(rank, lts+lte,'x')
    
    plot(rank,0*rank+mlta);
    errorbar(rank,lta,nerra,'o')
    
    lbs=['lin: ';'exp: ';'sum: '];
    ls(end+1)=mlta;
    le(end+1)=0;
    
    if gtexist
        [lsgt,ltsgt]=evaluatetime(nstart_gt,p,'m');
        [legt,ltegt]=evaluatetime(nend_gt,p,'g');
        lsgt(end+1)=0;
        legt(end+1)=0;
    
          plot(timepoints, ltsgt+ltegt,'d')
          
        title(axtt,[lbs num2str([ls le lsgt legt]*100,'%4.0f')]);
            legend('data start',['lin fit: ' num2str(ls(1)*100,'%2.0f')],...
       ['exp fit: ' num2str(ls(2)*100,'%2.0f')],'data end',...
       ['lin fit: ' num2str(le(1)*100,'%2.0f')],['exp fit: ' num2str(le(2)*100,'%2.0f')],...
       'sum',...
       'gt start',['lin fit: ' num2str(lsgt(1)*100,'%2.0f')],['exp fit: ' num2str(lsgt(2)*100,'%2.0f')], ...
       'gt end',['lin fit: ' num2str(legt(1)*100,'%2.0f')],['exp fit: ' num2str(legt(2)*100,'%2.0f')])
    
    else
        title(axtt,[lbs num2str([ls le]*100,'%4.1f, ')]);
            legend('data start',['lin fit: ' num2str(ls(1)*100,'%2.1f')],...
       ['exp fit: ' num2str(ls(2)*100,'%2.1f')],'data end',...
       ['lin fit: ' num2str(le(1)*100,'%2.1f')],['exp fit: ' num2str(le(2)*100,'%2.0f')])
    
    end
    ylabel('labeling efficiency')
 
    %number of localizations
    %aware of bias for not segmenting npcs with low number of localizations
axle=obj.initaxis('localizations');
hold off
hlocs=histogram(axle,numlocs,0:1:max(numlocs));


ns=floor(mean(numlocs)-std(numlocs));
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

if gtexist %not from simulation
ax6=obj.initaxis('GT filtered');
    p.ploton=false;
    bs_gtf=bootstrp(20,@fitNPClabeling,numcorners(indgood),p);
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

    ax6b=obj.initaxis('GT all');
        p.ploton=false;
    bs_gtfa=bootstrp(20,@fitNPClabeling,numcornersa(indgood),p);
    berr_gtfa=std(bs_gtfa);
    p.ploton=true;  
    
    hold off
    hnc=hist(numcornersa(indgood),nb);
    bar(nb,hnc)
    hold on
    pf=fitNPClabeling(hnc,p);
    results(1).gtall=pf*100; 
    results(2).gtall=berr_gtfa*100;
    title(['GT all: ' num2str(100*pf,'%2.1f') '\pm' num2str(100*berr_gtfa,'%2.1f')])
    axis tight
    numcorners=numcorners(indgood);
    rd=(0:length(numcorners)-1)/length(numcorners)/2;
ax7=obj.initaxis('comparison');
    plot(numcorners+rd,numfoundint-numcorners,'+')
    hold on
    plot(numcorners+rd,numfoundrat-numcorners,'x')
    plot(numcorners+rd,numbercornerassigneddirect -numcorners,'d')
    hold off
    legend('integer','prob','assigned d','all')

    numfoundintm=mean(numfoundint-numcorners);
    numfoundratm=mean(numfoundrat-numcorners);
    numbercornerassinedm=mean(numbercornerassinged-numcorners);
numbercornerassineddm=mean(numbercornerassigneddirect-numcorners);
    title(['int: ' num2str(numfoundintm,2),', rat: ' num2str(numfoundratm,2), ' ,asssigned d: ' num2str(numbercornerassineddm,2)])
    
    

sp=obj.locData.files.file(filefile).info.simulationParameters;
results(3).gapint=sp.labeling_efficiency;
results(3).gapf=sp.blinks;
results(3).assigned=sp.photons;
results(3).all=sp.lifetime;
results(3).gtfilt=sp.background;

results(3).gtall=0;
end
results(2).all=-1;

results(2).file=0;
results(3).file=0;
results(1).file=p.filecheck*(p.filenumbers);
% dat=struct2table(results);
axp=ax0.Parent;
delete(ax0);
ht=uitable('Parent',axp);
struct2uitable(ht, results,'flip','%2.1f')
% A=cell2mat(ht.Data);
copytoexcel(results,'flip');
ht.ColumnWidth={50,40,30,40};
ht.ColumnName={};
out=[];

if p.copy2page
    sm=3;
    sn=4;
    f=figure;
    f.Renderer='painters';
    ht2=ht.copy;
    ht2.Parent=f;
    axtemp= subplot(sm,sn,[1 2]);
    ht2.Units='normalized';
    ht2.Position=axtemp.Position;
    delete(axtemp)
%     subplot(2,3,1,axpsf2)
    if exist('axpsf','var')
        axpsf2=axpsf.copy;
        axpsf2.Parent=f;
        subplot(sm,sn,3,axpsf2)
%     axis tight
    end
    axt=ax1.copy;
    axt.Parent=f;
    subplot(sm,sn,4,axt)
    axt=ax2.copy;
    axt.Parent=f;
    subplot(sm,sn,8,axt)
        axt=ax3.copy;
    axt.Parent=f;
    subplot(sm,sn,7,axt)
    
    axttime=axtt.copy;
    axttime.Parent=f;
    subplot(sm,sn,[5 6 9 10],axttime)   
    
    axt=axle.copy;
    axt.Parent=f;
    subplot(sm,sn,11,axt)   
    
%         axt=ax4.copy;
%     axt.Parent=f;
%     subplot(sm,sn,10,axt)
    
    if exist('ax6','var')
        axt=ax6.copy;
    axt.Parent=f;
    subplot(sm,sn,11,axt) 
    end
    if exist('ax6b','var')
        axt=ax6b.copy;
    axt.Parent=f;
    subplot(sm,sn,12,axt) 
%      ax7l=legend(axt,t7,'Location','northoutside');
%     legend(axt,'integer','prob','assigned','Location','northwestoutside');
    end
    fn=se.files(p.filenumbers(1)).name;
    t=textwrap({fn},50);
    th=text(axttime,0.15,0.05,t,'FontSize',9,'Interpreter','none');
    
    
end

end


function [out,pf,nerr,tp]=evaluatetime(n,p,col)
if nargin <3
    col='k';
end
    
p.ploton=false;
p.fitrange=[1 8];
nerr=ones(1,size(n,2))/100;
pf=zeros(1,size(n,2));
for k=1:size(n,2)
%     h=hist(n(:,k),nb);
    nh=n(:,k);
    if sum(nh>0)>5 %minium for fitting)
        pf(k)=fitNPClabeling(nh,p);
        if p.bootstrap
            bs_numfoundint=bootstrp(10,@fitNPClabeling,nh,p);
            nerr(k)=std(bs_numfoundint);
        end
    else
        pf(k)=0;
        nerr(k)=0;
    end
end
% tp=timepoints(1,:);
tp=linspace(0,1,size(n,2)); %rank
% plot(tp,pf,'*');
% hold on
errorbar(tp,pf,nerr,[col '*'])
% fitrange=pf>.05 &pf<0.95;
fitrange=pf<0.95;
% fitrange=true(size(pf));

if p.bootstrap
    delta=min(nerr(nerr>0));
    w=1./(nerr+delta);
    w=w(fitrange);
    w=w/max(w);
else
    w=[];
end

fr=fit(tp(fitrange)',pf(fitrange)','poly1','Weights',w');

if pf(end)-pf(1) < 0 
    le=fr(0);
    o=tp(end);
    g = fittype( @(a,b,c,x) 1-a-c*exp(-b*(o-x)));
    fx=fit(tp(fitrange)',pf(fitrange)',g,'StartPoint',[pf(1) 1/tp(end) pf(1)],'Weights',w');
    lex=fx(0);
    
else
    fx=fit(tp(fitrange)',pf(fitrange)','1-a-c*exp(-b*x)','StartPoint',[pf(end) 1/tp(end) pf(end)],'Weights',w');
    le=fr(tp(end));
    lex=fx(tp(end));
end

hold on
plot(tp,fr(tp),col);
plot(tp,fx(tp),[col '--']);
xlabel('rank');
out=[le;lex];
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

pard.bootstrap.object=struct('String','Bootstrap error bars for time analysis','Style','checkbox');
pard.bootstrap.position=[4,1];
pard.bootstrap.Width=2;

pard.plugininfo.type='ROI_Analyze';
  

pard.copy2page.object=struct('String','Copy to own page','Style','checkbox');
pard.copy2page.position=[6,1];
pard.copy2page.Width=2;

pard.plugininfo.type='ROI_Analyze';

end