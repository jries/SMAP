classdef Segment_cotracks<interfaces.DialogProcessor
    %Links molecules in consecutive frames for SPT analysis
    methods
        function obj=Segment_cotracks(varargin)        
            obj@interfaces.DialogProcessor(varargin{:}) ;
            obj.showresults=true;
        end
        function out=run(obj,p)
            
            out=runi(obj,p);
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
       
    end
end
function out=runi(obj,p)
out=[];

% identify tracks in ch0 (min locs)
% longest track in ch1 with overlapping time.
% filter: locs ch1, progressive, fraction overlap...
% add to ROI manager
% option: ch1, ch2, co-tracks.


%filters


minlenlocs=p.minlenlocs;
cotracklength=p.cotracklength; %minimum data poits associated between tracks
cotrackfraction = p.cotrackfraction; %minimum fraction of localizations associated
lennmstartend=p.lennmstartend; % minimum endpoint- startpoint
layers=find(obj.getPar('sr_layerson'));
[locs,indin]=obj.locData.getloc({'xnm','ynm','time','tid','filenumber','thi'},'layer',layers,'position','all','grouping','ungrouped','removefilter','filenumber');


SE=obj.locData.SE;
if isempty(SE.cells)
    cg=ROIManager.Segment.makeCellGrid;
    cg.attachLocData(obj.locData);
    cg.attachPar(obj.P);
    p=cg.getAllParameters;
    cg.run(p);
end
cellp=vertcat(SE.cells.pos);
for c=length(SE.cells):-1:1
    cellfn(c,1)=SE.cells(c).info.filenumber;
end
% fhh=mode(locs.filenumber);
fnums=unique(locs.filenumber);
for ff=1:length(fnums)
    fhh=fnums(ff);

    tids0=locs.tid(locs.thi==0 & locs.filenumber==fhh);
    numloc=histcounts(tids0,1:max(tids0)+1);
    id0=find(numloc>minlenlocs);
    for k=1:length(id0)
        idh=locs.tid==id0(k) & locs.filenumber==fhh;
        xh=locs.xnm(idh); yh=locs.ynm(idh);
        time=locs.time(idh);
        len=sqrt((xh(end)-xh(1)).^2+(yh(end)-yh(1)).^2);
        if len<lennmstartend
            continue
        end
        %find other color
        ind1p=find(locs.time>=time(1) & locs.time<=time(end) & locs.thi==1 & locs.filenumber==fhh);
        if isempty(ind1p)
            continue
        end
        id1=mode(locs.tid(ind1p));
        numl1=sum(locs.tid==id1);
        if numl1<cotracklength || numl1/length(xh)< cotrackfraction
            continue
        end
    
        currentsite=interfaces.SEsites;
        currentsite.pos=[mean(xh), mean(yh), 0];    
        thisf=~(cellfn==fhh);
        [~,cind]=min(sum((currentsite.pos(1:2)-cellp).^2,2)+thisf*1e9);
        currentsite.info.cell=SE.cells(cind).ID;
        currentsite.info.filenumber=fhh;
        currentsite.annotation.comments=num2str(time([1 end])/1000);
        SE.addSite(currentsite);
    end
end
SE.processors.preview.updateSitelist;
end
             
function pard=guidef(obj)
% pard.showt.object=struct('String','Segment:','Style','text');
% pard.showt.position=[1,1];
% pard.showt.Width=1;
% pard.showtraces.object=struct('String',{{'co-tracks','ch0','ch1'}},'Style','popupmenu');
% pard.showtraces.position=[1,2];
% pard.showtraces.Width=1.;

pard.minlenlocst.object=struct('String','Min length track (locs)','Style','text');
pard.minlenlocst.position=[2,1];
pard.minlenlocst.Width=1.5;

pard.minlenlocs.object=struct('String','5','Style','edit');
pard.minlenlocs.position=[2,2.5];
pard.minlenlocs.Width=.5;

pard.lennmstartendt.object=struct('String','start-end (nm) >','Style','text');
pard.lennmstartendt.position=[2,3.5];
pard.lennmstartendt.Width=1.;
pard.lennmstartend.object=struct('String','300','Style','edit');
pard.lennmstartend.position=[2,4.5];
pard.lennmstartend.Width=.5;

pard.cotracklenght.object=struct('String','co track length (locs) >','Style','text');
pard.cotracklenght.position=[3,1];
pard.cotracklenght.Width=1.5;
pard.cotracklength.object=struct('String','4','Style','edit');
pard.cotracklength.position=[3,2.5];
pard.cotracklength.Width=.5;

pard.cotrackfractiont.object=struct('String','co track fraction >','Style','text');
pard.cotrackfractiont.position=[3,3.5];
pard.cotrackfractiont.Width=1.5;
pard.cotrackfraction.object=struct('String','0','Style','edit');
pard.cotrackfraction.position=[3,4.5];
pard.cotrackfraction.Width=.5;






pard.plugininfo.description=sprintf('co-tracking analysis');
pard.plugininfo.type='ProcessorPlugin';
end

