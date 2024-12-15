classdef SimpleTracking<interfaces.DialogProcessor
    %Links molecules in consecutive frames for SPT analysis
    methods
        function obj=SimpleTracking(varargin)        
            obj@interfaces.DialogProcessor(varargin{:}) ;
            obj.inputParameters={'sr_pixrec','numberOfLayers','sr_pos','sr_size','layers','sr_layerson'};
            obj.showresults=true;
        end
        function out=run(obj,p)
            
            out=tracki(obj,p);
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end
function out=tracki(obj,p)
switch p.trackwhat.selection
    case 'all'
        [locs,indin]=obj.locData.getloc({'xnm','ynm','znm','locprecnm','locprecznm','phot','frame','track_id','channel'},'position','roi','grouping','ungrouped');
        out=tracksingle(obj,p,locs,indin);
    case 'layers'
        layers=find(obj.getPar('sr_layerson'));
        for k=1:length(layers)
            [locs,indin]=obj.locData.getloc({'xnm','ynm','znm','locprecnm','locprecznm','phot','frame','track_id'},'layer',layers(k),'position','roi','grouping','ungrouped');
            out=tracksingle(obj,p,locs,indin);
            p.overwritetracks=false;
        end
    case 'channels'
        channels=unique(obj.locData.loc.channel);
        [locs,indin]=obj.locData.getloc({'xnm','ynm','znm','locprecnm','locprecznm','phot','frame','track_id','channel'},'position','roi','grouping','ungrouped');
        for k=1:length(channels)
            chgood=locs.channel==channels(k);
            out=tracksingle(obj,p,locs,indin,chgood);
            p.overwritetracks=false;
        end
end

tracklength=zeros(size(obj.locData.loc.xnm),'like',obj.locData.loc.xnm);
for k=1:max(obj.locData.loc.track_id)
    indt=obj.locData.loc.track_id==k;
    tracklength(indt)=sum(indt);
end
obj.locData.setloc('track_length',tracklength);


obj.locData.regroup;
obj.setPar('locFields',fieldnames(obj.locData.loc))


end

function out=tracksingle(obj,p,locs,indin,indfilter)
% locs.xnm=-min(locs.xnm);
% locs.ynm=-min(locs.ynm);
if nargin>4 %filter transmitted
    findin=find(indin);
    indin(findin(~indfilter))=false;
else
    indfilter=true(1,length(locs.xnm));
end
if ~isempty(locs.znm)
    xyzt=horzcat(locs.xnm(indfilter),locs.ynm(indfilter),locs.znm(indfilter),locs.frame(indfilter));
    p.dim=3;
else
    p.dim=2;
xyzt=double(horzcat(locs.xnm(indfilter),locs.ynm(indfilter),locs.frame(indfilter)));
end
p.quiet=false;
tracks=mytrack(xyzt,p.maxdisp,p);
llocs.x=locs.xnm(indfilter);llocs.y=locs.ynm(indfilter);llocs.frame=locs.frame(indfilter);

ftracks=tracks(:,end-1);
[~,trackssort]=sort(ftracks);
ltracks.x=tracks(trackssort,1);
ltracks.y=tracks(trackssort,2);
ltracks.frame=tracks(trackssort,end-1);
trackids=tracks(trackssort,end);
[iAa,iBa,nA,nB,nseen]=matchlocsall(llocs,ltracks,0,0,1);

% id=tracks(iBa,end);
id=trackids(iBa);
if isempty(locs.track_id)||p.overwritetracks
    idall=zeros(size(indin),'single');
    % tracklength=zeros(size(indin),'single');
    minID=0;
else
    idall=obj.locData.loc.track_id; 
    % tracklength=obj.locData.loc.track_length;
    minID=max(idall);
end
findin=find(indin);
idall(findin(iAa))=id+minID;
obj.locData.setloc('track_id',idall);


out=[];
% p.mode=1;
% analyze_SPT(tracks,p);
end
             
function pard=guidef
pard.t1.object=struct('String','Simple tracking','Style','text');
pard.t1.position=[1,2];
pard.t1.Width=4;

pard.maxdispt.object=struct('String','Maximum displacement between frames (nm)','Style','text');
pard.maxdispt.position=[2,1];
pard.maxdispt.Width=1.5;
pard.maxdispt.TooltipString='an estimate of the maximum distance that a particle would move in a single time interval.';

pard.maxdisp.object=struct('String','1000','Style','edit');
pard.maxdisp.position=[2,2.5];
pard.maxdisp.Width=.5;
pard.maxdisp.TooltipString=pard.maxdispt.TooltipString;

pard.memt.object=struct('String','Maximum dark time (frames)','Style','text');
pard.memt.position=[2,3];
pard.memt.Width=1.5;
pard.memt.TooltipString=sprintf(['this is the number of time steps that a particle can be \n'...
            'lost and then recovered again.  If the particle reappears \n'...
           'after this number of frames has elapsed, it will be \n'...
            'tracked as a new particle. The default setting is zero. \n'...
           'this is useful if particles occasionally drop out of the data.']);
pard.mem.object=struct('String','2','Style','edit');
pard.mem.position=[2,4.5];
pard.mem.Width=.5;
pard.mem.TooltipString=pard.memt.TooltipString;


pard.goodt.object=struct('String','Minimum length of tracks','Style','text');
pard.goodt.position=[3,1];
pard.goodt.Width=1.5;
pard.goodt.TooltipString=sprintf(['set this keyword to eliminate all trajectories with \n'...
            ' fewer than param.good valid positions.  This is useful \n'...
           'due to blinking noise particles in the data stream.']);
pard.good.object=struct('String','2','Style','edit');
pard.good.position=[3,2.5];
pard.good.Width=.5;
pard.good.TooltipString=pard.goodt.TooltipString;

pard.trackwhat.object=struct('String',{{'all','layers','channels'}},'Style','popupmenu');
pard.trackwhat.position=[3,3];
pard.trackwhat.Width=1.5;
% pard.uselayers.object=struct('String','use layers','Style','checkbox','Value',0);
% pard.uselayers.position=[4,1];
% pard.uselayers.Width=1.5;

% pard.usechannels.object=struct('String','use channels','Style','checkbox','Value',1);
% pard.usechannels.position=[4,3];
% pard.usechannels.Width=1.5;

pard.overwritetracks.object=struct('String','overwrite tracks','Style','checkbox','Value',1);
pard.overwritetracks.position=[5,1];
pard.overwritetracks.Width=1.5;
% pard.overwritetracks.TooltipString=pard.goodt.TooltipString;


pard.plugininfo.description=sprintf('Links molecules in consecutive frames for SPT analysis');
pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.Name='uTrack';
end

