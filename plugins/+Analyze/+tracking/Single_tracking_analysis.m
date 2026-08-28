classdef Single_tracking_analysis<interfaces.DialogProcessor
    %Links molecules in consecutive frames for SPT analysis
    properties
        xyaxis
        tracks
    end
    methods
        function obj=Single_tracking_analysis(varargin)        
            obj@interfaces.DialogProcessor(varargin{:}) ;
            obj.showresults=true;
        end
        function out=run(obj,p)
            
            out=runi(obj,p);
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function loadtiff(obj,a,b,field,sstr)
            oldf=obj.getSingleGuiParameter(field);
            if isempty(oldf)
                oldf=sstr;
            end

            [newf,newpath]=uigetfile(oldf);
            if newf
                obj.setGuiParameters(struct(field,[newpath newf]));
            end
        
        end
    end
end
function out=runi(obj,p)
out=[];
% set filter values

minlenframes=p.minlenframes;


% test directed movement
aspectratio=p.aspectratio; %should be below this for directed motion.
lennmstartend=p.lennmstartend; % minimum endpoint- startpoint

minvelocity=p.velocitymin; %nm/s
maxvelocity=p.velocitymax;

%analyze everything displayed
layers=find(obj.getPar('sr_layerson'));
obj.locData.filter; %does this fix the bug?
[locs,indin]=obj.locData.getloc({'xnm','ynm','znm','xpix','ypix','frame','track_id','track_length','layer','filenumber','channel'},'layer',layers,'position','roi','grouping','ungrouped');
if isempty(locs.track_id)
    error('please use the SimpleTracking plugin first')
end

exposuretime=obj.locData.files.file(locs.filenumber(1)).info.exposure;
roi=obj.locData.files.file(locs.filenumber(1)).info.roi;
pixelsize=obj.locData.files.file(locs.filenumber(1)).info.cam_pixelsize_um*1000;

% coarse analysis
% test directed movement, calculate statistics
usetracks=unique(locs.track_id(locs.track_id>0));
trackstat.lennmstartend=zeros(max(usetracks),1);
trackstat.stdshort=zeros(max(usetracks),1);
trackstat.angle=zeros(max(usetracks),1);
trackstat.stdlong=zeros(max(usetracks),1);
trackstat.lenframe=zeros(max(usetracks),1);
trackstat.velocity=zeros(max(usetracks),1);
locs.track_length_new=zeros(size(locs.xnm));
trackstat.processive=false(max(usetracks),1);
trackstat.runtime=zeros(max(usetracks),1);
trackstat.runlength=zeros(max(usetracks),1);
trackstat.x0=zeros(max(usetracks),1);
trackstat.y0=zeros(max(usetracks),1);
trackstat.channel=zeros(max(usetracks),1);

for k=1:length(usetracks)
    iduset=usetracks(k);
    tind=find(locs.track_id==iduset);
    trackstat.lenframe(iduset)=length(tind);
    locs.track_length_new(tind)=trackstat.lenframe(iduset);
    if trackstat.lenframe(iduset)<2
        continue
    end
    
    %length
    xh=locs.xnm(tind);yh=locs.ynm(tind);
    trackstat.lennmstartend(iduset)=sqrt((xh(end)-xh(1)).^2+(yh(end)-yh(1)).^2);
    
    %aspect ratio
    c = cov(xh-mean(xh), yh-mean(yh));
    [a, ev] = eig(c);
    [ev,ind] = sort(diag(ev));
    [xa, ya] = deal(a(1,ind(end)), a(2,ind(end)));
    trackstat.angle(iduset)=cart2pol(xa, ya);
    trackstat.stdshort(iduset)=real(sqrt(ev(1)));
    trackstat.stdlong(iduset)=real(sqrt(ev(2)));
    trackstat.channel(iduset)=mode(locs.channel(tind));
end
intrack=(locs.track_id>0);
longtracks=locs.track_length_new>minlenframes;

usetracks=unique(locs.track_id(intrack&longtracks)); %these are the interesting clusters.

%look for processive movement
validstats= trackstat.lenframe>minlenframes;
validstats = validstats & trackstat.stdshort./trackstat.stdlong<aspectratio;
validstats = validstats & trackstat.lennmstartend > lennmstartend;
trackstat.filter=validstats;

% plot tracks
ax=obj.initaxis('xy');
obj.xyaxis=ax;
cols=[1 1 0
      0 1 1];
colgood=[0 0 0];
msize=3;
lw=1;
symb='-';
colind=(trackstat.filter)+1;   

if p.plottracks
for k=1:length(usetracks)  
    idh=usetracks(k);
    indtr=locs.track_id==idh;
    hp=plot(ax,locs.xnm(indtr)/pixelsize(1)-roi(1),locs.ynm(indtr)/pixelsize(2)-roi(2),symb,'Color',cols(colind(idh),:),'LineWidth',lw,'Tag','test','MarkerSize',msize);
    dtRows = [dataTipTextRow("frame",double(locs.frame(indtr))),...
        dataTipTextRow("ID",double(locs.track_id(indtr)))];
        alldatatip=vertcat(hp.DataTipTemplate.DataTipRows,dtRows');   
        hp.DataTipTemplate.DataTipRows=alldatatip;
    hold(ax,"on")
end
axis(ax,'equal');
xlim(ax,[0 roi(3)])
ylim(ax,[0 roi(4)/2])
axis(ax,'ij');
end

% Analyze processivity
processiveids=find(trackstat.filter);
if contains(p.showtraces.selection,'processive') 
    figure; 
    f=0;
end
v=zeros(length(processiveids),1);runlength=v; runtime=v; rmse=v; gof=v;channel=v;
goodv=true(length(processiveids),1);
plotind=1;
% statvall=struct;
for k=1:length(processiveids)
    idt=locs.track_id==processiveids(k);
     % idt=locs.track_id==285;
    statv=analyse_processive_track(locs.xnm(idt),locs.ynm(idt),locs.frame(idt));
    tmin=(min(locs.frame(idt)));
    tmax=(max(locs.frame(idt)));
    exposuretimes=exposuretime/1000;
    % filter by analysis results
    if statv.gof>p.gofmax || statv.rmse > p.rmsemax || statv.numpoints<minlenframes || statv.v/exposuretimes<minvelocity
        goodv(k)=false;
    end
    
    v(k)=statv.v/exposuretimes;
    statv.vnms=v(k);
    trackstat.velocity(processiveids(k))=v(k);
    rmse(k)=statv.rmse; gof(k)=statv.gof;
    runlength(k)=statv.len;
    trackstat.runlength(processiveids(k))=runlength(k);
    runtime(k)=statv.numpoints;
    channel(k)=trackstat.channel(processiveids(k));
    trackstat.runtime(processiveids(k))=runtime(k);
    trackstat.processive(processiveids(k))=goodv(k);
    trackstat.x0(processiveids(k))=statv.x0;
    trackstat.y0(processiveids(k))=statv.y0;
    % trackstat.t0(processiveids(k))=tmin; %no, rather tmin
    ff='%2.0f';
    %plot good tracks
    if contains(p.showtraces.selection,'processive') && goodv(k)
        if 2*plotind-f>30
            f=f+30;
            figure
        end
        ax1=subplot(5,6,2*plotind-1-f);
        ax2=subplot(5,6,2*plotind-f);    
        plotind=plotind+1;
        plot(ax1,statv.plot.time,statv.plot.xr,statv.plot.time,statv.plot.xfitr);
        xlabel(ax1,'time')
        ylabel(ax1,'xrot')
        plot(ax2,statv.plot.x,statv.plot.y,statv.plot.xfit,statv.plot.yfit);
        axis(ax2,'equal')
        xlabel(ax2,'xrot')
        ylabel(ax2,'yrot')
        title(ax2,['f: ' num2str(tmin) ':', num2str(tmax),', x: ' num2str(mean(locs.xpix(idt)),'%3.0f') ', y: ' num2str(mean(locs.ypix(idt)),'%3.0f')])   
        title(ax1,"Id:"+processiveids(k)+ ", v: " + num2str(v(k),ff) + ", fit: " + num2str(gof(k),'%3.1f')+ "," + num2str(rmse(k),'%3.0f'))
    end    

    %plot tracks to overview
    if goodv(k) & p.plottracks
        hp=plot(ax,locs.xnm(idt)/pixelsize(1)-roi(1),locs.ynm(idt)/pixelsize(2)-roi(2),symb,'Color',colgood,'LineWidth',lw,'Tag','test','MarkerSize',msize);
        hold(ax,"on")
        dtRows = [dataTipTextRow("frame",double(locs.frame(idt))),...
        dataTipTextRow("ID",double(locs.track_id(idt)))];
        alldatatip=vertcat(hp.DataTipTemplate.DataTipRows,dtRows');   
        hp.DataTipTemplate.DataTipRows=alldatatip;
    end
    statvall(processiveids(k))=statv;
end
if ~exist('statvall','var')
    statvall=[];
end
% trackstat.processive=processiveids;
out.tracks=trackstat;
out.processive=trackstat.processive;
out.processiveids=processiveids;
out.trackids=usetracks;
obj.tracks.trackstat=trackstat;
obj.tracks.analysis=statvall;


axv=obj.initaxis('v');
histogram(axv,v(goodv))
title(axv, "v nm/s, mean "+mean(v(goodv)) + ", median " + median(v(goodv)));
axr=obj.initaxis('runlength');
histogram(axr,runlength(goodv))
title(axr,"length nm, mean "+mean(runlength(goodv)) + ", median " + median(runlength(goodv)));
axrt=obj.initaxis('runtime');
histogram(axrt,runtime(goodv))
title(axrt,"length time frames, mean "+mean(runtime(goodv)) + ", median " + median(runtime(goodv)));

axv=obj.initaxis('rmse');
hold(axv,'off')
h1=histogram(axv,rmse); 
hold(axv,'on');
h2=histogram(axv,rmse(goodv));
h1.BinWidth=20;h2.BinWidth=20;
axv=obj.initaxis('gof');
h1=histogram(axv,gof); 
hold(axv,'on');
h2=histogram(axv,gof(goodv));
h1.BinWidth=.2;h2.BinWidth=.2;

% extracts filename from file path:
fnum=find(obj.getPar('sr_layerson'),1,'first');
fl=obj.getPar('filelist_long').String;
filePath=string(fl{fnum});

% filePath = string(obj.getPar('lastSMLFile'));
tablename=strrep(filePath,'_sml.mat','_tracks.csv'); %if this does not work, replace fileslm by filePath
[path,fileName]=fileparts(filePath); 

output=(table(string(fileName),sum(trackstat.lenframe>=minlenframes), sum(goodv), mean(v(goodv)), mean(runlength(goodv)), mean(runtime(goodv)),'VariableNames', {'Filename','Total', 'processive','Velocity','runlength','runtime'}));

disp(output)
out.summarytable=output;

outputtracks=(table(repmat(string(fileName),sum(goodv),1), processiveids(goodv),(v(goodv)), (runlength(goodv)), (runtime(goodv)), channel(goodv),'VariableNames', {'Filename','ID','Velocity','runlength','runtime','channel'}));
% writetable(outputtracks,tablename) % creates a CSV output for each input
out.trackstable=outputtracks;

end

function plotcurrent_callback(a,b,obj)
f=obj.xyaxis.Parent.Parent.Parent;
dcm=datacursormode(f);
c=getCursorInfo(dcm);
if isempty(c)
    disp('select a track with datacursor')
    return
end
ID=c(1).Target.DataTipTemplate.DataTipRows(end).Value(1);
% x=c.Target.XData; y=c.Target.YData;
% angle=obj.tracks.trackstat.angle(ID);
% [xr,yr]=rotcoord(x-mean(x),y-mean(y),angle);

fit=obj.tracks.analysis(ID);
plt=fit.plot;
figure(788)
subplot(1,2,1);
plot(plt.time,plt.xr,'b',fit.filtered.time,fit.filtered.xr,'b+',plt.time,plt.xfitr,'r')

title("v " + num2str(fit.vnms,'%2.0f')+ ", len " + num2str(fit.numpoints))
subplot(1,2,2);
plot(plt.xr,plt.yr,'b',fit.filtered.xr,fit.filtered.yr,'b+',plt.xfitr,plt.yfitr,'r')
axis equal
title("gof " + num2str(fit.gof,'%2.1f')+ ", rmse " + num2str(fit.rmse,'%2.0f'))

end
             
function pard=guidef(obj)

pard.minlenframest.object=struct('String','Min length track (frames)','Style','text');
pard.minlenframest.position=[1,1];
pard.minlenframest.Width=1.5;

pard.minlenframes.object=struct('String','10','Style','edit');
pard.minlenframes.position=[1,2.5];
pard.minlenframes.Width=.5;

pard.dmt.object=struct('String','directed movement:','Style','text');
pard.dmt.position=[2,1];

pard.aspectratiot.object=struct('String','aspect ratio <','Style','text');
pard.aspectratiot.position=[2,2];
pard.aspectratiot.Width=1.;
pard.aspectratio.object=struct('String','1','Style','edit');
pard.aspectratio.position=[2,2.7];
pard.aspectratio.Width=.5;


pard.lennmstartendt.object=struct('String','start-end (nm) >','Style','text');
pard.lennmstartendt.position=[2,3.5];
pard.lennmstartendt.Width=1.;
pard.lennmstartend.object=struct('String','200','Style','edit');
pard.lennmstartend.position=[2,4.5];
pard.lennmstartend.Width=.5;


pard.velocityt.object=struct('String','velocity (nm/s) min','Style','text');
pard.velocityt.position=[4,1];
pard.velocityt.Width=1.5;
pard.velocitymin.object=struct('String','200','Style','edit');
pard.velocitymin.position=[4,2.5];
pard.velocitymin.Width=.5;
pard.velocitymt.object=struct('String','max','Style','text');
pard.velocitymt.position=[4,4];
pard.velocitymt.Width=0.5;
pard.velocitymax.object=struct('String','1500','Style','edit');
pard.velocitymax.position=[4,4.5];
pard.velocitymax.Width=.5;


pard.goft.object=struct('String','Goodness of Fit rmse/diff(rmse) <','Style','text');
pard.goft.position=[5,1];
pard.goft.Width=2;
pard.gofmax.object=struct('String','5','Style','edit');
pard.gofmax.position=[5,3];
pard.gofmax.Width=.5;
pard.rmset.object=struct('String','rmse (nm) <','Style','text');
pard.rmset.position=[5,3.5];
pard.rmset.Width=1.;
pard.rmsemax.object=struct('String','100','Style','edit');
pard.rmsemax.position=[5,4.5];
pard.rmsemax.Width=.5;

pard.showt.object=struct('String','Show tracks:','Style','text');
pard.showt.position=[7,1];
pard.showt.Width=0.5;
pard.showtraces.object=struct('String',{{'none','processive'}},'Style','popupmenu');
pard.showtraces.position=[7,1.5];
pard.showtraces.Width=1.5;

pard.plottracks.object=struct('String','plot tracks','Style','checkbox','Value',1);
pard.plottracks.position=[7,3];
pard.plottracks.Width=1;

pard.plotcurrent.object=struct('String','plot selected','Style','pushbutton','Callback',{{@plotcurrent_callback,obj}});
pard.plotcurrent.position=[7,4];
pard.plotcurrent.Width=1;

pard.plugininfo.description=sprintf('co-tracking analysis');
pard.plugininfo.type='ProcessorPlugin';
end

