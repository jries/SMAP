classdef DCoffset_Beads_ROIs<interfaces.DialogProcessor&interfaces.SEProcessor
    properties
    end
    methods
        function obj=DCoffset_Beads_ROIs(varargin)        
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
out=[];
p.mainch=0;
p.followch=~p.mainch;
switch p.mode.selection
    case 'fast'
        fastoffset(obj,p);   
    case 'tails'
        tails2c(obj,p);
    case 'drift'
        driftoffset(obj,p);
    case 'FoV'
        driftfov(obj,p);
end
end

function fastoffset(obj,p)
p.axx=obj.initaxis('x');
p.axy=obj.initaxis('y');
p.axdx=obj.initaxis('dx');
p.axdy=obj.initaxis('dy');
sites=obj.SE.sites;
layers=find(obj.getPar('sr_layerson'));
subx=p.subtractx*p.subtractmed;suby=p.subtracty*p.subtractmed;

for k=1:length(sites)
    if ~sites(k).annotation.use
        continue
    end
    % mainch=0;
    followch=~p.mainch;
    locs=obj.locData.getloc({'xnm','ynm','time','tid','thi'},'layer',layers,'Position',sites(k),'grouping','ungrouped');
    indmain=locs.thi==p.mainch;
    tid1=mode(locs.tid(indmain));
    ind1=locs.tid==tid1;
    id2=findch2(locs, ind1, followch);
    ind2=locs.tid==id2;

    x1h=locs.xnm(ind1)-subx;x2h=locs.xnm(ind2);y1h=locs.ynm(ind1)-suby;y2h=locs.ynm(ind2);
    t1h=locs.time(ind1);t2h=locs.time(ind2);

    dx=distancedc(t1h, x1h, t2h, x2h);
    dy=distancedc(t1h, y1h, t2h, y2h);

   
    
    plot(p.axdx, t1h(p.skip+1:end), dx(p.skip+1:end))
    plot(p.axx,t1h(p.skip+1:end), x1h(p.skip+1:end), t2h(p.skip+1:end), x2h(p.skip+1:end))
    
    plot(p.axdy, t1h(p.skip+1:end), dy(p.skip+1:end))
    plot(p.axy, t1h(p.skip+1:end), y1h(p.skip+1:end), t2h(p.skip+1:end), y2h(p.skip+1:end))

    hold(p.axdx,"on")
    hold(p.axdy,"on")
    hold(p.axx,"on")
    hold(p.axy,"on")

    xm(k)=mean(dx); ym(k)=mean(dx);
    stx(k)=std(dx);sty(k)=std(dy);
   
end
xlabel(p.axx,'time (ms)'); ylabel(p.axx,'x (nm)')
xlabel(p.axy,'time (ms)'); ylabel(p.axy,'y (nm)')
xlabel(p.axdx,'time (ms)'); ylabel(p.axdx,'dx (nm)')
xlabel(p.axdy,'time (ms)'); ylabel(p.axdy,'dy (nm)')

title(p.axdx,['dx: ' num2str(xm,'%2.1f, ') ' std: ' num2str(stx, '%2.1f, ')] )
title(p.axdy,['dy: ' num2str(ym,'%2.1f, ') ' std: ' num2str(sty, '%2.1f, ')] )
end

function driftoffset(obj,p)
p.axdx=obj.initaxis('dx');
p.axdy=obj.initaxis('dy');
p.axdxn=obj.initaxis('dx norm');
p.axdyn=obj.initaxis('dy norm');
p.axarrow=obj.initaxis('distance');

p.t1=obj.initaxis('dt start vs stop');
p.t2=obj.initaxis('dt fraction vs start');


% p.dyims=obj.initaxis('dy s');
colorsarrow=jet(256)*0.9;
maxt=max(obj.locData.loc.time);
subx=p.subtractx*p.subtractmed;suby=p.subtracty*p.subtractmed;
layers=find(obj.getPar('sr_layerson'));

sites=obj.SE.sites;
for ss=length(sites):-1:1
    if ~sites(ss).annotation.use
        continue
    end
    followch=~p.mainch;
    locs=obj.locData.getloc({'xnm','ynm','time','tid','thi','ecc','eco'},'layer',layers,'Position',sites(ss),'grouping','ungrouped');
    indmain=locs.thi==p.mainch;
    if isempty(locs.tid)
        continue
    end
    hc=histcounts(locs.tid(indmain),1:(max(locs.tid))+1);
    
    tidgood=find(hc>=p.minlen);
    if isempty(tidgood)
        continue
    end
    % hc(tidgood)

    dx=zeros(length(tidgood),1);dy=dx; ti=dx;x1=dx;x2=dx; y1=dx;y2=dy;
    dtstart=dx; dtstop=dx; tfraction=dx;
 
    for k=length(tidgood):-1:1
        ind1=locs.tid==tidgood(k);
        id2=findch2(locs, ind1, followch);
        ind2=locs.tid==id2;
        if sum(ind2)==0
            dx(k)=NaN; dy(k)=NaN;
            dtstart=NaN; dtstop=NaN; tfraction=NaN;
            continue
        end
        x1h=locs.xnm(ind1)-subx;x2h=locs.xnm(ind2);y1h=locs.ynm(ind1)-suby;y2h=locs.ynm(ind2);
        t1h=locs.time(ind1);t2h=locs.time(ind2);
        % mx1=median(x1h(p.skip+1:end));
        % mx2=median(x2h(p.skip+1:end));
        % my1=median(y1h(p.skip+1:end));
        % my2=median(y2h(p.skip+1:end));
        dx(k)=median(distancedc(t1h(p.skip+1:end), x1h(p.skip+1:end), t2h(p.skip+1:end), x2h(p.skip+1:end)));
        dy(k)=median(distancedc(t1h(p.skip+1:end), y1h(p.skip+1:end), t2h(p.skip+1:end), y2h(p.skip+1:end)));
        ti(k)=t1h(p.skip);
        x1(k)=median(x1h(p.skip+1:end)); x2(k)=median(x2h(p.skip+1:end));
        y1(k)=median(y1h(p.skip+1:end)); y2(k)=median(y2h(p.skip+1:end));
        x1a(ss,k)=x1(k);y1a(ss,k)=y1(k);dxa(ss,k)=dx(k);dya(ss,k)=dy(k);
        trel=ceil(t1h(1)/maxt*255);
        plot(p.axarrow, [x1(k) x1(k)+dx(k)* p.vectorlength],[y1(k) y1(k)+dy(k)* p.vectorlength],'Color',colorsarrow(trel,:));
        hold(p.axarrow,"on")

        dtstart(k)=t2h(1)-t1h(1);
        dtstarta(ss,k)=dtstart(k);
        dtstop(k)=t1h(end)-t2h(end);
        tfraction(k)=(t2h(end)-t2h(1))/(t1h(end)-t1h(1));
    end
    linest=p.plotline;
        plot(p.axdx, ti, dx,linest); hold(p.axdx,"on")
        plot(p.axdy, ti, dy, linest); hold(p.axdy,"on")
        plot(p.axdxn, ti, dx-dx(1),linest,ti,0*ti,'k'); hold(p.axdxn,"on")
        plot(p.axdyn, ti, dy-dy(1), linest,ti,0*ti,'k'); hold(p.axdyn,"on")
        plot(p.t1, dtstart+(rand(length(dtstart),1))*0.25,dtstop+(rand(length(dtstart),1))*0.25, linest); hold(p.t1,"on")
        plot(p.t2, dtstart+(rand(length(dtstart),1))*0.25,tfraction, linest); hold(p.t2,"on")
end

% scatter(x1a(:),y1a(:),5,dxa(:),'filled','Parent',p.dxims); colorbar(p.dxims)
% scatter(x1a(:),y1a(:),5,dya(:),'filled','Parent',p.dyims); colorbar(p.dyims)
% plot(p.dxims,dxa(:))

xlabel(p.axdxn,'time (ms)'); ylabel(p.axdxn,'dx-dx(1) (nm)')
xlabel(p.axdyn,'time (ms)'); ylabel(p.axdyn,'dy-dy(1) (nm)')
xlabel(p.axdx,'time (ms)'); ylabel(p.axdx,'dx (nm)')
title(p.axdx,"median: "+ median(abs(dxa(dx~=0)))+", std: "+num2str(std(dxa(dx~=0)),'%2.1f') + " nm")
xlabel(p.axdy,'time (ms)'); ylabel(p.axdy,'dy (nm)')
% title(p.axdy,median(abs(dya(dx~=0))))
title(p.axdy,"median: "+ median(abs(dya(dx~=0)))+", std: "+num2str(std(dya(dx~=0)),'%2.1f') + " nm")
xlabel(p.t1,'dt start (ms)'); ylabel(p.t1,'dt stop (ms)')
xlabel(p.t2,'dt start (ms)'); ylabel(p.t2,'fraction ch1/ch0')

xlim(p.t1,quantile(dtstarta(dtstarta>0),[0. .98]))
xlim(p.t2,quantile(dtstarta(dtstarta>0),[0. .98]))

colormap(p.axarrow,'jet')
colorbar(p.axarrow)
title(p.axarrow,'d (t)')
dxx=quantile(dxa(dx~=0),[.99]);
xlim(p.axarrow,quantile(x1a(x1a~=0),[.01 .99])+[-1 1]*dxx*p.vectorlength)
ylim(p.axarrow,quantile(y1a(y1a~=0),[.01 .99])+[-1 1]*dxx*p.vectorlength)
end


function driftfov(obj,p)
p.dxim=obj.initaxis('dx');
p.dyim=obj.initaxis('dy');
p.dxims=obj.initaxis('dx s');
p.dyims=obj.initaxis('dy s');
% p.axdxn=obj.initaxis('dx norm');
% p.axdyn=obj.initaxis('dy norm');

sites=obj.SE.sites;
dx=zeros(length(sites),1);dy=dx; x1=dx;x2=dx;y1=dx; y2=dx; ti=dx;
layers=find(obj.getPar('sr_layerson'));
subx=p.subtractx*p.subtractmed;suby=p.subtracty*p.subtractmed;
for k=1:length(sites)
    if ~sites(k).annotation.use
        continue
    end
    followch=~p.mainch;
    locs=obj.locData.getloc({'xnm','ynm','time','tid','thi'},'layer',layers,'Position',sites(k),'grouping','ungrouped');
    indmain=locs.thi==p.mainch;
    tid1=mode(locs.tid(indmain));
    ind1=find(locs.tid==tid1);
    id2=findch2(locs, ind1, followch);
    ind2=find(locs.tid==id2);
    if length(ind1)<p.minlen || length(ind2)<p.minlen
        continue
    end

    x1h=locs.xnm(ind1)-subx;x2h=locs.xnm(ind2);y1h=locs.ynm(ind1)-suby;y2h=locs.ynm(ind2);
    t1h=locs.time(ind1);t2h=locs.time(ind2);
    dx(k)=median(distancedc(t1h(p.skip+1:end), x1h(p.skip+1:end), t2h(p.skip+1:end), x2h(p.skip+1:end)));
    dy(k)=median(distancedc(t1h(p.skip+1:end), y1h(p.skip+1:end), t2h(p.skip+1:end), y2h(p.skip+1:end)));
    ti(k)=t1h(p.skip);
    x1(k)=median(x1h(p.skip+1:end)); x2(k)=median(x2h(p.skip+1:end));
    y1(k)=median(y1h(p.skip+1:end)); y2(k)=median(y2h(p.skip+1:end));

end
[X,Y]=meshgrid(min(x1):100:max(x1),min(y1):100:max(y1));
dxs=scatteredInterpolant(x1,y1,dx);
imagesc(p.dxim, dxs(X,Y));
colorbar(p.dxim)
dys=scatteredInterpolant(x1,y1,dy);
imagesc(p.dyim, dys(X,Y));
colorbar(p.dyim)

scatter(x1,y1,50,dx,'filled','Parent',p.dxims); colorbar(p.dxims)
scatter(x1,y1,50,dy,'filled','Parent',p.dyims); colorbar(p.dyims)
end

function tails2c(obj,p)
p.meantails=obj.initaxis('meantails');
p.tails1=obj.initaxis('tails ch0');
p.tails2=obj.initaxis('tails ch1');

p.taillength1=obj.initaxis('taillength ch0');
p.taillength2=obj.initaxis('taillength ch1');
p.taildurations=obj.initaxis('tailduration ch0 - ch1');
layers=find(obj.getPar('sr_layerson'));
sites=obj.SE.sites;
tailsall1=[];tailind=1; phot1=[];tailsall2=[]; phot2=[];

for ss=1:length(sites)
    if ~sites(ss).annotation.use
        continue
    end
    followch=~p.mainch;
    locs=obj.locData.getloc({'xnm','ynm','time','tid','thi','ecc','eco'},'layer',layers,'Position',sites(ss),'grouping','ungrouped');
    indmain=locs.thi==p.mainch;
    hc=histcounts(locs.tid(indmain),1:max(locs.tid)+1);
    tidgood=find(hc>=p.minlen);

    dx=zeros(length(tidgood),1);dy=dx; ti=dx;x1=dx;x2=dx; y1=dx;y2=dy;taillen1=dx;taillen2=dx;taillocs1=dx; taillocs2=dx;
    for k=1:length(tidgood)
        ind1=locs.tid==tidgood(k);
        id2=findch2(locs, ind1, followch);
        ind2=locs.tid==id2;
        if sum(ind2)==0
            dx(k)=NaN; dy(k)=NaN;
            continue
        end
        x1h=locs.xnm(ind1);x2h=locs.xnm(ind2);y1h=locs.ynm(ind1);y2h=locs.ynm(ind2);
        t1h=locs.time(ind1);t2h=locs.time(ind2);
        mx1=median(x1h(p.skip+1:end));
        mx2=median(x2h(p.skip+1:end));
        my1=median(y1h(p.skip+1:end));
        my2=median(y2h(p.skip+1:end));
        dx(k)=median(distancedc(t1h(p.skip+1:end), x1h(p.skip+1:end), t2h(p.skip+1:end), x2h(p.skip+1:end)));
        dy(k)=median(distancedc(t1h(p.skip+1:end), y1h(p.skip+1:end), t2h(p.skip+1:end), y2h(p.skip+1:end)));
        ti(k)=t1h(p.skip);
        x1(k)=median(x1h(p.skip+1:end)); x2(k)=median(x2h(p.skip+1:end));
        y1(k)=median(y1h(p.skip+1:end)); y2(k)=median(y2h(p.skip+1:end));

        d1=sqrt((locs.xnm(ind1)-x1(k)).^2 +(locs.ynm(ind1)-y1(k)).^2);
        taillen1(k)=max(d1);
        d2=sqrt((locs.xnm(ind2)-x2(k)).^2 +(locs.ynm(ind2)-y2(k)).^2);
        taillen2(k)=max(d2);
        % if tailind<p.maxplottails
        %     plot(p.tails,  sqrt((x1h-mx1).^2+(y1h-my1).^2));hold(p.tails,"on")
        % end

        % tl1=find(d1>p.prec,1,'last');
        tl1=find(d1<p.prec,1,'first');
        if ~isempty(tl1)
            taillocs1(k)=tl1;
        end
       
        % tl2=find(d2>p.prec,1,'last');
        tl2=find(d2<p.prec,1,'first');

        if ~isempty(tl2)
            taillocs2(k)=tl2;
        end
        tailsall1(1:length(d1),tailind)=d1;
        phot1(1:length(d1),tailind)=locs.ecc(ind1)+locs.eco(ind1);
        tailsall2(1:length(d2),tailind)=d2;
        phot2(1:length(d2),tailind)=locs.ecc(ind2)+locs.eco(ind2);           
        tailind=tailind+1;
    end
    plot(p.taillength1,taillocs1+ (rand(length(taillocs1),1))*0.5,taillen1,'.');  hold(p.taillength1,"on")
    plot(p.taillength2,taillocs2+ (rand(length(taillocs2),1))*0.5,taillen2,'.');  hold(p.taillength2,"on")
    plot(p.taildurations,taillocs1+ (rand(length(taillocs1),1))*0.5,taillocs2+ (rand(length(taillocs2),1))*0.5,'.');  hold(p.taildurations,"on")
   
end
tailsall1(tailsall1==0)=NaN;tailsall2(tailsall2==0)=NaN;
phot1(phot1==0)=NaN;phot2(phot2==0)=NaN;
phottot1=cumsum(mean(phot1,2,'omitnan'));
phottot2=cumsum(mean(phot2,2,'omitnan'));
plot(p.meantails,phottot1,mean(tailsall1,2,'omitnan'),'-+')
hold(p.meantails,"on")
plot(p.meantails,phottot2,mean(tailsall2,2,'omitnan'),'-o')
plot(p.meantails, phottot1, 120./sqrt(phottot1),'k')
legend(p.meantails,'MINFLUX Ch0','MINFLUX Ch1','PSF/sqrt(N)')
xlabel(p.meantails,'total photons'); ylabel(p.meantails,'localization error (nm)')
p.meantails.YAxis.Scale="log";

xlabel(p.taillength1,'duration of tails (localizations)'); ylabel(p.taillength1,'length of tails (nm)')
xlabel(p.taillength2,'duration of tails (localizations)'); ylabel(p.taillength2,'length of tails (nm)')
plotevery=size(tailsall1,2)/p.maxplottails;
plot(p.tails1, tailsall1(:,1:plotevery:end));
xlabel(p.tails1,'localization'); ylabel(p.tails1,'distance to average final position (nm)')
plot(p.tails2, tailsall2(:,1:plotevery:end));
xlabel(p.tails2,'localization'); ylabel(p.tails2,'distance to average final position (nm)')


xlabel(p.taildurations,'duration of tails ch0 (localizations)'); ylabel(p.taildurations,'duration of tails ch1 (localizations)')
end




function dx=distancedc(t1, x1, t2,x2)
if isempty(x2)
    dx=NaN;
    return;
end
dx=0*t1;
indt2=1;
    for indt1=1:length(t1)
        while indt2<length(t2) && t2(indt2)<t1(indt1)
            indt2=indt2+1;
        end
        if indt2>1 && t2(indt2)-t1(indt1)>t1(indt1)-t2(indt2-1)
            indt2=indt2-1;
        end
        dx(indt1)=x1(indt1)-x2(indt2);
    end
end

     
function id2=findch2(locs, ind1, followch, maxd, offset)
tm=min(locs.time(ind1)); tx=max(locs.time(ind1));
id2=mode(locs.tid(locs.time > tm & locs.time < tx & locs.thi==followch));
end

function pard=guidef(obj)
pard.modet.object=struct('String','analysis mode','Style','text');
pard.modet.position=[1,1];
pard.mode.object=struct('String',{{'fast','drift','FoV','tails'}},'Style','popupmenu');
pard.mode.position=[1,2];

pard.subtractmed.object=struct('String','subtract pos [x,y] (nm)','Style','checkbox','Value',true);
pard.subtractmed.position=[2,1];
pard.subtractmed.Width = 2;
pard.subtractx.object=struct('String','0','Style','edit');
pard.subtractx.position=[2,3];
pard.subtractx.Width=0.5;
pard.subtracty.object=struct('String','0','Style','edit');
pard.subtracty.position=[2,3.5];
pard.subtracty.Width=0.5;


pard.skipt.object=struct('String','skip first','Style','text');
pard.skipt.position=[3,1];
pard.skip.object=struct('String','20','Style','edit');
pard.skip.position=[3,2];

pard.minlt.object=struct('String','min length','Style','text');
pard.minlt.position=[3,3];
pard.minlen.object=struct('String','50','Style','edit');
pard.minlen.position=[3,4];

pard.prect.object=struct('String','precision (nm)','Style','text');
pard.prect.position=[4,1];
pard.prec.object=struct('String','4','Style','edit');
pard.prec.position=[4,2];

pard.maxplottailst.object=struct('String','max number of plots','Style','text');
pard.maxplottailst.position=[4,3];
pard.maxplottails.object=struct('String','100','Style','edit');
pard.maxplottails.position=[4,4];

pard.vectorlengtht.object=struct('String','factor distance arrow plot','Style','text');
pard.vectorlengtht.position=[5,1];
pard.vectorlength.object=struct('String','10','Style','edit');
pard.vectorlength.position=[5,2];

pard.plotlinet.object=struct('String','Line style','Style','text');
pard.plotlinet.position=[5,3];
pard.plotline.object=struct('String','-','Style','edit');
pard.plotline.position=[5,4];

pard.plugininfo.type='ROI_Analyze';
pard.plugininfo.description=' ';

end