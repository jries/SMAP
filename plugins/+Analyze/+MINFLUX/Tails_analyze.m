classdef Tails_analyze<interfaces.DialogProcessor&interfaces.SEProcessor
    properties
    end
    methods
        function obj=Tails_analyze(varargin)        
            obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={};
            obj.showresults=true;
        end
        
        function out=run(obj,p)  
            tic
            out=runintern(obj,p);
            toc
            
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end

    end
end


function out=runintern(obj,p)
out=[];

p.meantails=obj.initaxis('meantails');
p.meantailsadd=obj.initaxis('meantailsadd','keep');
p.tails=obj.initaxis('tails');
p.taillengthduration=obj.initaxis('length vs duration');
p.taillength=obj.initaxis('length');
p.duration=obj.initaxis('duration');

layers=find(obj.getPar('sr_layerson'));

tidgood=[];
maxtid=max(obj.locData.loc.tid);
if contains(p.source.selection, 'ROI') %use ROI manager
% make a collection of tids
    for ss=1:length(obj.SE.sites)
        if ~obj.SE.sites(ss).annotation.use
            continue
        end
        locs=obj.locData.getloc({'xnm','ynm','time','tid','filenumber'},'layer',layers,'Position',obj.SE.sites(ss),'grouping','ungrouped');
        tidfn=locs.tid+maxtid*double(locs.filenumber);
        hc=histcounts(tidfn,1:max(tidfn)+1);
        tidgoodh=find(hc>=p.minlen);
        tidgood(end+1:end+length(tidgoodh))=tidgoodh;
    end
else %look at all valid IDs
    locs=obj.locData.getloc({'xnm','ynm','time','tid','filenumber'},'layer',layers,'Position','all','grouping','ungrouped');
    tidfn=locs.tid+maxtid*double(locs.filenumber);
    hc=histcounts(tidfn,1:max(tidfn)+1);
    tidgood=find(hc>=p.minlen);
end

tailsall=[];phot=[];
tailsall(1,length(tidgood))=0;phot(1,length(tidgood))=0;

tskip=zeros(length(tidgood),1);xskipmed=zeros(length(tidgood),1); yskipmed=zeros(length(tidgood),1);
taillen=zeros(length(tidgood),1);taillocs=zeros(length(tidgood),1);

locs=obj.locData.getloc({'xnm','ynm','time','tid','ecc','eco','filenumber'},'layer',layers,'Position','all','grouping','ungrouped');
tidfn=locs.tid+maxtid*double(locs.filenumber);
locsall=obj.locData.getloc({'xnm','ynm','tid','ecc','eco','time','filenumber'},'Position','all','grouping','ungrouped');
tidfnall=locsall.tid+maxtid*double(locsall.filenumber);

if contains(p.source.selection, 'ROI')
locs=locsall; tidfn=tidfnall;
end
for k=1:length(tidgood)
    ind=find(tidfn==tidgood(k));

    xh=locs.xnm(ind);yh=locs.ynm(ind);
    th=locs.time(ind);

    tskip(k)=th(p.skip);
    endl=min(length(xh),p.maxlen);
    xskipmed(k)=median(xh(p.skip+1:endl)); 
    yskipmed(k)=median(yh(p.skip+1:endl)); 

    dist=sqrt((xh(1:endl)-xskipmed(k)).^2 +(yh(1:endl)-yskipmed(k)).^2);
    taillen(k)=max(dist);

    tl1=find(dist<p.prec,1,'first');
    % tl1=find(dist>p.prec,1,'last');
    if ~isempty(tl1)
        taillocs(k)=tl1;
    else
        taillocs(k)=NaN;
    end
    tailsall(1:length(dist),k)=dist;
    
    %
    inda=find(tidfnall==tidgood(k));
    ino=1;
    for m=1:endl
        phota=0;
        while ino<=length(inda) && locsall.time(inda(ino))<=th(m)
            phota = phota + locsall.ecc(inda(ino))+locsall.eco(inda(ino));
            ino=ino+1;
        end
        phot(m,k)=phota;
    end
    % photo(1:length(dist),k)=locs.ecc(ind)+locs.eco(ind);
end


% plot(p.taillengthduration,taillocs+ (rand(length(taillocs),1))*0.5,taillen,'.');  hold(p.taillengthduration,"on")
axes(p.taillengthduration); hold(p.taillengthduration,'on');
plot(p.taillengthduration, taillocs + rand(size(taillocs))*0.5, taillen, '.');


ylim([0, quantile(taillen,.99)])

n=0:5:quantile(taillen,.99);
histogram(taillen,n,'Parent',p.taillength)
xlabel(p.taillength,'length of tail (nm)')
hTL = histogram(taillen,n,'Parent',p.taillength); legend(p.taillength, hTL, ['median = ', num2str(round(median(taillen, "omitnan"),1))]);
n = 1:max(taillocs, [], "omitnan") + 1;
histogram(taillocs,n,'Parent',p.duration)
xlabel(p.duration,'localizations until convergence')
% plot(p.taillength,h)
hTL = histogram(taillocs,n,'Parent',p.duration); legend(p.duration, hTL, ['median = ', num2str(median(taillocs, "omitnan"))]);


tailsall(tailsall==0)=NaN;
phottot=cumsum(mean(phot,2,'omitnan'));


plot(p.meantails,phottot,mean(tailsall,2,'omitnan'),'-+')
hold(p.meantails,"on")
plot(p.meantails, phottot, 120./sqrt(phottot),'k')
legend(p.meantails,'MINFLUX','PSF/sqrt(N)')
xlabel(p.meantails,'total photons'); ylabel(p.meantails,'localization error (nm)')
p.meantails.YAxis.Scale="log";

xlabel(p.taillengthduration,'duration of tails (localizations)'); ylabel(p.taillengthduration,'length of tails (nm)')

plotevery=length(tidgood)/p.maxplottails;
hold(p.tails,'off')
plot(p.tails, tailsall(:,1:plotevery:end));
hold(p.tails,'on')
plot(p.tails, mean(tailsall,2),'k','LineWidth',2);
xlabel(p.tails,'localization'); ylabel(p.tails,'distance to average final position (nm)')

[~,ftagname]=fileparts(strrep(obj.locData.files.file(locs.filenumber(1)).name,'\',filesep));

hold(p.meantailsadd,'on')
plot(p.meantailsadd, mean(tailsall,2),'DisplayName',ftagname);
xlabel(p.meantailsadd,'localization'); ylabel(p.meantailsadd,'distance to average final position (nm)')
legend(p.meantailsadd);
end


function pard=guidef(obj)
pard.sourcet.object=struct('String','Localizations from: ','Style','text');
pard.sourcet.position=[1,1];
pard.source.object=struct('String',{{'All filtered','ROI manager'}},'Style','popupmenu');
pard.source.position=[1,2];

pard.skipt.object=struct('String','skip first','Style','text');
pard.skipt.position=[3,1];
pard.skip.object=struct('String','20','Style','edit');
pard.skip.position=[3,2];

pard.minlt.object=struct('String','min length','Style','text');
pard.minlt.position=[3,3];
pard.minlen.object=struct('String','50','Style','edit');
pard.minlen.position=[3,4];



pard.prect.object=struct('String','Convergence: ','Style','text');
pard.prect.position=[4,1];
pard.precmode.object=struct('String',{{'first d < ', 'last d > '}},'Style','popupmenu');
pard.precmode.position=[4,2];
pard.prec.object=struct('String','2','Style','edit');
pard.prec.position=[4,3];
pard.prec.Width=0.5;
pard.prect2.object=struct('String','nm','Style','text');
pard.prect2.position=[4,3.5];



pard.maxplottailst.object=struct('String','max number of plots','Style','text');
pard.maxplottailst.position=[5,1];
pard.maxplottails.object=struct('String','100','Style','edit');
pard.maxplottails.position=[5,2];

pard.maxlt.object=struct('String','max length analysis','Style','text');
pard.maxlt.position=[5,3];
pard.maxlen.object=struct('String','50','Style','edit');
pard.maxlen.position=[5,4];

pard.plugininfo.type='ROI_Analyze';
pard.plugininfo.description=' ';

end