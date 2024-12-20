classdef StepsMINFLUX<interfaces.SEEvaluationProcessor
%     
    properties
        peakline
        roihandle
        axstep
        axxy
        xyplotoffset
        steps
        coord
        range
        index
        
        stats
        id;
        locsuse;
        manualcuration
        steps_x_range
        xy_range
    end
    methods
        function obj=StepsMINFLUX(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function out=run(obj, p)
           obj.steps=[];obj.range=[];obj.coord=[];
           obj.steps_x_range=[];
           obj.xy_range=[];
           if isfield(obj.site.evaluation,obj.name) 
                out=obj.site.evaluation.(obj.name);
                if isfield(obj.site.evaluation.(obj.name),'steps')
                    obj.steps=out.steps;
                end
                if isfield(obj.site.evaluation.(obj.name),'range')
                    obj.range=out.range;
                else
                    obj.range=[];
                end                
           end

           %identify all localizations in track
           usefields={'xnm','ynm','groupindex','tid','time','znm', 'efo', 'cfr', 'eco', 'ecc', 'efc','filenumber','vld','sta'};
           locs=obj.getLocs(usefields,'layer',find(obj.getPar('sr_layerson')),'size',obj.getPar('se_siteroi')/2,'removeFilter',{'time'});
           if isempty(locs.xnm)
                disp('no localizations')
                return
           end
           if contains(p.link.selection,'group')
               fid='groupindex';
           else
               fid='tid';
           end
            
           times=str2num(obj.site.annotation.comments)*1000;
           if isempty(times)
               id=mode(locs.(fid));
           else
               timemin=obj.site.image.parameters.layerparameters{1}.colorfield_min;
               timemax=obj.site.image.parameters.layerparameters{1}.colorfield_max;
               if max(times)<1000
                   times=timemin+(timemax-timemin)*times/1000;
               end
               
               if length(times)==1 % one element
                   timewinselect=(timemax-timemin)*0.05;
                   times=times+[-timewinselect timewinselect];
               end
               ind=locs.time > times(1) & locs.time < times(2);
               id=mode(locs.(fid)(ind));
               
           end

           if isfield(p,'fromdc') && p.fromdc
               switch p.dcch
                   case 1
                       idname='id1';   
                   case 2
                       idname='id2';
               end
               if ~isfield(obj.site.evaluation.StepsMINFLUX_dc,idname)
                   out=[];
                   return
               end
               id=obj.site.evaluation.StepsMINFLUX_dc.(idname);
               % obj.locData.SE.sites(1).evaluation
           end
           filenumberh=mode(locs.filenumber);
           
           if p.filterlocs
               locs=obj.locData.getloc(usefields,'layer',find(obj.getPar('sr_layerson')),'position','all','removeFilter',{'time','filenumber'});
               obj.locsuse=locs;
               index=obj.locsuse.(fid) == id & obj.locsuse.filenumber ==filenumberh;
           else
               obj.locsuse=obj.locData.loc;
               index=obj.locsuse.(fid)==id & obj.locsuse.filenumber ==filenumberh;
           end
           if p.onlyvld
                index=index & obj.locsuse.vld==1;
           end
           obj.index=index;
           obj.id=id;
           
            
           indexvld=index & obj.locsuse.vld==1;
     
           % XX decide if to get filtered or real coordinates 
           %get coordinates
           x=obj.locsuse.xnm(indexvld);
           y=obj.locsuse.ynm(indexvld);
           time=obj.locsuse.time(indexvld);

           if isfield(obj.locData.loc,'znm')
               z=obj.locsuse.znm(indexvld);
           else
               z=[];
           end

           [xr,yr,angle]=rotateCenterCoordinates(x,y,time,obj.range);
            % if p.filtertrackmode
                zf=z;
                windowsize=p.filterwindow;
                xf=runningWindowAnalysis(time,x,time,windowsize,p.filtermode.selection);     
                yf=runningWindowAnalysis(time,y,time,windowsize,p.filtermode.selection);  
                if ~isempty(z)
                    zf=runningWindowAnalysis(time,z,time,windowsize,p.filtermode.selection);  
                end
%                 plot(axx,timeplot,xf,'b');
                [xfr,yfr,angle]=rotateCenterCoordinates(xf,yf,time,obj.range,angle);
            % else
                
                % xfr=[]; yfr=[]; zf=[];
            % end


           

%            c = cov(x-mean(x), y-mean(y));
%            [a, ev] = eig(c);
%            [ev,ind] = sort(diag(ev));
%            [xa, ya] = deal(a(1,ind(end)), a(2,ind(end)));
%            angle = cart2pol(xa, ya);
%            [xr,yr]=rotcoord(x-mean(x),y-mean(y),angle);
%            indx=xr<mean(xr);
%            if mean(time(indx))>mean(time(~indx)) %increasing position with time
%                xr=-xr;
%            end
            
%            xr=x;yr=y;  %XXXX to not rotate

           obj.coord.xr=xr;obj.coord.yr=yr;obj.coord.time=time;obj.coord.timeplot=time-min(time);
           obj.coord.x=x;obj.coord.y=y;obj.coord.xfr=xfr;obj.coord.yfr=yfr;
           obj.coord.z=z;  obj.coord.zf=zf;
          
           if  isempty(obj.steps) || p.refitalways
               refit(0,0,obj,1)
           end
           calculatestepparameters(obj, obj.steps.indstep);
           plotsteps(obj)
           out=obj.site.evaluation.(obj.name);

           out.statall=calculatestatistics(obj,index);
           out.stattrack=calculatestatistics(obj,index,obj.coord.indtime);

           filelist=obj.getPar('filelist_short');
           out.filename=filelist.String{mode(obj.locsuse.filenumber)};
           
           plotstatistics(obj)
        end
        function updateroiposition(obj,a,b)
            %catch delete and add vertex!!!
            oldtime=b.PreviousPosition(:,1);
            newtime=b.CurrentPosition(:,1);
           movedvertex= find(oldtime ~= newtime);
            timemoved=newtime(movedvertex);
            [~,newind]=min((obj.coord.timeplot(obj.coord.indtime)-timemoved).^2);
            if (movedvertex== 1 || newind>obj.steps.indstep(movedvertex-1)) && (movedvertex== length(obj.steps.indstep) ||newind<obj.steps.indstep(movedvertex+1))
                obj.steps.indstep(movedvertex)=newind;
                obj.steps.steptime(movedvertex)=timemoved;
                calculatestepparameters(obj,obj.steps.indstep)
            end
            manualcurate(obj)
        end
        function addvertexroi(obj,a,b)
            %catch delete and add vertex!!!
            oldtime=obj.steps.steptime;
            newtime=a.Position(:,1);
            ind=find(oldtime~=newtime(1:end-1),1,'first');
            if isempty(ind)
                [~,newind]=min((obj.coord.timeplot(obj.coord.indtime)-newtime(end)).^2);
                newsteps = [obj.steps.indstep;newind];
            else
                [~,newind]=min((obj.coord.timeplot(obj.coord.indtime)-newtime(ind)).^2);
                newsteps=[obj.steps.indstep(1:ind-1); newind ;obj.steps.indstep(ind:end)];
            end
            calculatestepparameters(obj,newsteps)
            manualcurate(obj)
        end
        function deletevertexroi(obj,a,b)
            oldtime=obj.steps.steptime;
            newtime=a.Position(:,1);
            ind=find(oldtime(1:end-1)~=newtime,1,'first');
            if isempty(ind)
                newsteps=obj.steps.indstep(1:end-1);
            else
                newsteps=[obj.steps.indstep(1:ind-1); obj.steps.indstep(ind+1:end)];
            end

            calculatestepparameters(obj,newsteps)
            manualcurate(obj)    
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end     
    end

end


function plotstatistics(obj)
index=obj.index;
% time=obj.locData.loc.time(index);
time=obj.locsuse.time(index);
dt=diff(time);
dtmin=min(dt(dt>0));
dtmedian=median(dt);
dtmean=mean(dt);
ltime=time-min(time);
ff='%2.1f';
axt=obj.setoutput('time');
bins=dtmin/2:dtmin:quantile(dt,0.995);
if length(bins)<3
    bins=dtmin/2:dtmin/2:dtmin*2;
end
hold(axt,'off')
histogram(axt,dt,bins)
hold(axt,'on')
histogram(axt,dt,bins)
xlabel(axt,'dt (ms)')
ylabel(axt,'frequency')
title(axt,['dtmin = ' num2str(dtmin,ff) ' ms, dtmedian = ' num2str(dtmedian,ff) ' ms, dtmean = ' num2str(dtmean,ff) ' ms.'])

axdt=obj.setoutput('dt');
hold(axdt,'off')
plot(axdt,ltime(2:end),dt)
if ~isempty(obj.range)
    hold(axdt,'on')
    yr=[min(dt) max(dt)];
    plot(axdt,obj.range(1)*[1 1],yr,'k--','HitTest','off')
    plot(axdt,obj.range(2)*[1 1],yr,'k--','HitTest','off')
end

xlabel(axdt,'time (ms)')
ylabel(axdt,'dt (ms)')


plotsimple(obj,'efo')
plotsimple(obj,'cfr')
plotsimple(obj,'eco')
plotsimple(obj,'ecc')
plotsimple(obj,'sta')

if obj.getSingleGuiParameter('msdanalysis') %MSD
    amsd=obj.setoutput('MSD');
    hold(amsd,'off')
    msd=obj.stats.msd;
    plot(amsd,msd.dt,msd.mean,'b',msd.dt,msd.mean+msd.std,'c',msd.dt,msd.mean-msd.std,'c')
    xlabel(amsd,'dt(ms)')
    ylabel(amsd,'msd (nm^2)')
    hold(amsd,'on')
    plot(amsd,msd.tfit,msd.fit(msd.tfit),'r');
    title(amsd,['D = ' num2str(msd.D) ' nm^2/ms, offset = ' num2str(msd.off) ' nm']);
end
end

function plotsimple(obj,field)
if ~isempty(obj.locsuse.(field))
    time=obj.locsuse.time(obj.index);
    ltime=time-min(time);
    axe=obj.setoutput(field);
    yval=obj.locsuse.(field)(obj.index);
    hold(axe,'off')
    plot(axe,ltime,yval)
    hold(axe,'on')
    if ~isempty(obj.range)
        yr=[min(yval) max(yval)];
        plot(axe,obj.range(1)*[1 1],yr,'k--','HitTest','off')
        plot(axe,obj.range(2)*[1 1],yr,'k--','HitTest','off')
    end
    xlabel(axe,'time (ms)')
    ylabel(axe,field)
    vall=obj.site.evaluation.(obj.name).statall.(field);
    vt=obj.site.evaluation.(obj.name).stattrack.(field);
    title(axe,[field ': ' num2str(vall,'%2.3g') ', ' num2str(vt,'%2.3g') ' (track)'])

end
end

function out=calculatestatistics(obj,indexin,indind)
if islogical(indexin)
    indexin=find(indexin);
end
if nargin<=2
    indind=true(size(indexin));
end
index=indexin(indind);
time=obj.locsuse.time(index);

dt=diff(time);
out.dtmin=min(dt);
out.dtmedian=median(dt);
out.dtmean=mean(dt);

out.id=obj.id;

out.efo=median(obj.locsuse.efo(index));
out.cfr=median(obj.locsuse.cfr(index));
out.eco=median(obj.locsuse.eco(index));
out.ecc=median(obj.locsuse.ecc(index));
out.efc=median(obj.locsuse.efc(index));
out.sta=median(obj.locsuse.sta(index));
out.nlocs=length((index));
out.tracktime=max(time)-min(time);
% xh=obj.coord.xr(indind);
xh=obj.coord.xr;
out.tracklength=xh(end)-xh(1);
out.velocity=out.tracklength/out.tracktime;

out.stdall.x=sqrt(sum(obj.steps.std.x.^2.*obj.steps.numlocsstep)/sum(obj.steps.numlocsstep));
out.stdall.y=sqrt(sum(obj.steps.std.y.^2.*obj.steps.numlocsstep)/sum(obj.steps.numlocsstep));
out.stdalldet.x=sqrt(sum(obj.steps.stddetrend.x.^2.*obj.steps.numlocsstep)/sum(obj.steps.numlocsstep));
out.stdalldet.y=sqrt(sum(obj.steps.stddetrend.y.^2.*obj.steps.numlocsstep)/sum(obj.steps.numlocsstep));

out.stdalldetmean.x=(sum(obj.steps.stddetrend.x.*obj.steps.numlocsstep)/sum(obj.steps.numlocsstep));
out.stdalldetmean.y=(sum(obj.steps.stddetrend.y.*obj.steps.numlocsstep)/sum(obj.steps.numlocsstep));

% if ~isempty(obj.coord.xfr)
    out.stdall.xf=sqrt(sum(obj.steps.std.xf.^2.*obj.steps.numlocsstep)/sum(obj.steps.numlocsstep));
    out.stdall.yf=sqrt(sum(obj.steps.std.yf.^2.*obj.steps.numlocsstep)/sum(obj.steps.numlocsstep));
    out.stdalldet.xf=sqrt(sum(obj.steps.stddetrend.xf.^2.*obj.steps.numlocsstep)/sum(obj.steps.numlocsstep));
    out.stdalldet.yf=sqrt(sum(obj.steps.stddetrend.yf.^2.*obj.steps.numlocsstep)/sum(obj.steps.numlocsstep));
    out.stdalldetmean.xf=(sum(obj.steps.stddetrend.xf.*obj.steps.numlocsstep)/sum(obj.steps.numlocsstep));
    out.stdalldetmean.yf=(sum(obj.steps.stddetrend.yf.*obj.steps.numlocsstep)/sum(obj.steps.numlocsstep));
% end


if obj.getSingleGuiParameter('msdanalysis') %MSD
    x=obj.coord.xr; %here no smoothing used. add option for this as well?
    y=obj.coord.yr;
    t=obj.coord.timeplot;
    
    msdc=(x-x').^2+(y-y').^2;
    dt=abs(t-t');
    tres=2;
    tmax=100;
    tmaxplot=min(max(dt(:),tmax*5));
    tb=(0:tres:tmaxplot)';
    [dtsort,indsort]=sort(dt(:));
    msd.mean=bindata(dtsort,msdc(indsort),tb,'mean');
    msd.std=bindata(dtsort,msdc(indsort),tb,'std');
    msd.dt=tb;
    
    indmax=round(tmax/tres);
    msd.mfit=msd.mean(1:indmax);
    msd.tfit=tb(1:indmax);
    msd.fit=fit(double(msd.tfit),double(msd.mfit),'poly1');
    msd.D=msd.fit.p1/4;
    msd.off=sqrt(msd.fit.p2);
    % hold(amsd,'on')
    % plot(amsd,tfit,fmsd(tfit));
    % title(amsd,['D = ' num2str(fmsd.p1/4) ' nm^2/ms, offset = ' num2str(sqrt(fmsd.p2)) ' nm']);
    out.msd=msd;
end

obj.stats=out;
end

function refit(a,b,obj,what)
p=obj.getAllParameters;
%stepfinder
if  ~isempty(obj.range)
    indtime=obj.coord.timeplot>=obj.range(1) & obj.coord.timeplot<=obj.range(2);
else
    indtime=true(size(obj.coord.xr));
end
obj.coord.indtime=indtime;

if strcmp(obj.getSingleGuiParameter('filtertrackmode').selection,'smooth stepfind')
   xr=obj.coord.xr(indtime);
else
   xr=obj.coord.xfr(indtime);
end

% p.fitmean=contains(p.fitmode.selection,'mean');
if what==1
    p.stepfunction=p.fitmode.selection;
    indstep=findstepsMINFLUX(xr,p);
%     try
%     steps=AutoStepfinderRies(obj.coord.xr(indtime),p);
%     catch err
%         steps=[];
%     end
%     if isempty(steps)
%         disp('no steps found')
%         steps.properties.StepSize=10;
%         steps.properties.IndexStep=[];
%     end
%     indstep=[1 ;steps.properties.IndexStep ];      
%     
%     if obj.getSingleGuiParameter('splitmerge')
%         stepsize=p.splitmergestep;
%         if isempty(stepsize)
%             stepsize=median(steps.properties.StepSize);
%         end
%         [indstep,stepvalue]=splitmergefit(obj.coord.xr(indtime),stepsize,p,steps);
%     end
else %refine
     [svalfit, indstep]=fitstepind(xr,obj.steps.indstep,str2func(p.fitmode.selection));
end

calculatestepparameters(obj,indstep);

plotsteps(obj)
end

function plotsteps(obj)
ff='%2.1f';
try
dcm_obj = datacursormode(obj.axstep.Parent.Parent.Parent);
info=dcm_obj.getCursorInfo;
catch err
    info=[];
end
ax2=obj.setoutput('steps_x');
hold(ax2,'off')

if ~isempty(obj.steps_x_range)
    obj.steps_x_range=ax2.XLim;
end

trackid=mode(obj.locsuse.tid(obj.index));

plotfiltered=obj.getSingleGuiParameter('filtertrackmode').Value>1 ;
if plotfiltered
    colornotfiltered=[1 1 1]*0.6;
    plot(ax2,obj.coord.timeplot,obj.coord.xr,'LineWidth',0.25,'Color',colornotfiltered,'HitTest','off');
    hold(ax2,'on')
    plot(ax2,obj.coord.timeplot,obj.coord.xfr,'k','HitTest','off');
    
    sxdetrend=std(diff(obj.coord.xr))/sqrt(2);sxfdetrend=std(diff(obj.coord.xfr))/sqrt(2);
    title(ax2,['std(x) detrend = ' num2str(sxdetrend,ff) ' nm.' ' std(xf) detrend = ' num2str(sxfdetrend,ff) ' nm.' ' std_s(x) = ' num2str(obj.stats.stdall.x,ff) ' nm.' ' std_sf(x) = ' num2str(obj.stats.stdall.xf,ff) ' nm.'])

else
    colornotfiltered='k';
    plot(ax2,obj.coord.timeplot,obj.coord.xr,'k','HitTest','off');
    hold(ax2,'on')
end


xlabel(ax2,'time (ms)')
ylabel(ax2,'position (nm)')
if ~isempty(obj.range)
    plot(ax2,obj.range(1)*[1 1],[min(obj.coord.xr) max(obj.coord.xr)],'k--','HitTest','off')
    plot(ax2,obj.range(2)*[1 1],[min(obj.coord.xr) max(obj.coord.xr)],'k--','HitTest','off')
    tm=obj.range(2);
else
    tm=max(obj.coord.timeplot);
end

legend(ax2,['tid: ' num2str(trackid)])

mv=obj.steps.stepvalue;
hstep=stairs(ax2,[obj.steps.steptime ;tm],[mv ;mv(end)],'r');
cm = uicontextmenu(ax2.Parent.Parent.Parent);
m = uimenu(cm,'Text','Menu1');
cm.ContextMenuOpeningFcn = @(src,event)disp('Context menu opened');
hstep.ContextMenu = cm;
dmv=diff(mv);
% for k=1:length(dmv)
%     text(ax2,(obj.steps.steptime(k+1)),mean(mv(k:k+1)),num2str(dmv(k),'%2.1f'),'FontSize',10,'Color','magenta','HitTest','off')
% end

showvalues(obj,ax2,obj.steps.steptime(2:end),(mv(1:end-1)+mv(2:end))/2,dmv);

if isempty(obj.steps_x_range)
    obj.steps_x_range=ax2.XLim;
else
    ax2.XLim=obj.steps_x_range;
end
obj.axstep=ax2;

if ~isempty(info)
    datatip(hstep,info.Position(1),info.Position(2));
end


%xy plot
goff=median(mod(obj.steps.stepvalue,16),'omitnan');
obj.xyplotoffset=goff;
% goff=0; %switch off shift 
axxy=obj.setoutput('xy');

if ~isempty(obj.xy_range)
    obj.xy_range=axxy.XLim;
end

hold(axxy,'off')
plot(axxy, obj.coord.xr-goff, obj.coord.yr,'Color',colornotfiltered)
axis(axxy,'equal')
xlabel(axxy,'x (nm)')
ylabel(axxy,'y (nm)')
hold(axxy,'on')

if plotfiltered
    plot(axxy, obj.coord.xfr-goff, obj.coord.yfr,'k')
end
obj.axxy=axxy;

scatter(axxy,obj.steps.stepvalue-goff,obj.steps.possteps.y,'r')
grid(axxy,'on')
axm=-16:-16:axxy.XLim(1);
axxy.XTick=[axm(end:-1:1) 0:16:axxy.XLim(2)];
axxy.YTick=round((axxy.YLim(1):6:axxy.YLim(2))/6)*6;

if isempty(obj.xy_range)
    obj.xy_range=axxy.XLim;
else
    axxy.XLim=obj.xy_range;
end

%plot selection box
if ~isempty(obj.range)
indmin=find(obj.coord.timeplot>=obj.range(1),1,'first');
indmax=find(obj.coord.timeplot<=obj.range(2),1,'last');
ylims=obj.axxy.YLim;
plot(axxy,(obj.coord.xr(indmin)-goff)*[1 1],ylims,'k--')
plot(axxy,(obj.coord.xr(indmax)-goff)*[1 1],ylims,'k--')
end

sigmax=std(obj.coord.xr);sigmay=std(obj.coord.yr);
sxdetrend=std(diff(obj.coord.xr))/sqrt(2);sydetrend=std(diff(obj.coord.yr))/sqrt(2);
[~, sxrobust]=robustMean(obj.coord.xr); [~, syrobust]=robustMean(obj.coord.yr);
title(axxy,['std(x) = ' num2str(sigmax,ff) ' nm, std(x) detrend = ' num2str(sxdetrend,ff) ' nm.' ' std(y) = ' num2str(sigmay,ff) ' nm, std(y) detrend = ' num2str(sydetrend,ff) ' nm.'])



if ~isempty(obj.coord.z)
    ax3D=obj.setoutput('zxy');
    hold(ax3D,'off')
    plot3(ax3D,obj.coord.xr-goff, obj.coord.yr,obj.coord.z-median(obj.coord.z))
    xlabel(ax3D,'x (nm)')
    ylabel(ax3D,'y (nm)')
    zlabel(ax3D,'z (nm)')
    hold(ax3D,'on')
    scatter3(ax3D,obj.steps.stepvalue-goff,obj.steps.possteps.y,obj.steps.possteps.z-median(obj.coord.z),'k')
    view(ax3D,0,0)
    axis(ax3D,'equal')
    
    axm=-16:-16:ax3D.XLim(1);
    ax3D.XTick=[axm(end:-1:1) 0:16:ax3D.XLim(2)];
%     ax3D.YTick=round((ax3D.YLim(1):25:ax3D.YLim(2))/25)*25;
%     ax3D.ZTick=round((ax3D.ZLim(1):25:ax3D.ZLim(2))/25)*25;
        ax3D.YTick=round((min(obj.coord.yr):25:max(obj.coord.yr))/25)*25;
    ax3D.ZTick=round((min(obj.coord.z-median(obj.coord.z)):25:max(obj.coord.z-median(obj.coord.z)))/25)*25;
%     ax3D.ZLim=[min(obj.coord.z-median(obj.coord.z)) max(obj.coord.z-median(obj.coord.z))];
%     ax3D.DataAspectRatio=[1 1 1];

    grid(ax3D,'on')
    szdetrend=std(diff(obj.coord.z))/sqrt(2);
    title(ax3D,['std_z(detrend) = ' num2str(szdetrend,ff)]);
end

axsy=obj.setoutput('steps_y');hold(axsy,'off')
plot(axsy,obj.coord.timeplot,obj.coord.yr,'Color',colornotfiltered)
if plotfiltered
    hold(axsy,'on')
    plot(axsy,obj.coord.timeplot,obj.coord.yfr,'k')
end

axsy.YTick=round((axxy.YLim(1):6:axxy.YLim(2))/6)*6;
try
axsy.XTick=obj.steps.steptime;
catch
end
axsy.XTickLabel=round(obj.steps.steptime);

axsy.GridColor='r';
grid(axsy);

if ~isempty(dmv) && length(dmv)>1
ax3=obj.setoutput('stephist');
histogram(ax3,dmv,min(dmv):.5:max(dmv));
%step time
axstept=obj.setoutput('dwelltime');
dt=max(0.1,round(mean(obj.steps.dwelltime)/5));
histogram(axstept,obj.steps.dwelltime,0:dt:max(obj.steps.dwelltime))
title(axstept,"mean step time = "+ num2str(mean(obj.steps.dwelltime),'%2.1f') + " ms");
end

%correlation
axcc=obj.setoutput('correlation');
x=obj.coord.xr(obj.coord.indtime);
h=histcounts(x,min(x):1:max(x));
xc=myxcorr(h,h);
plot(axcc,xc);
xlabel(axcc,'delta x (nm)');
ylabel(axcc,'auto corr')

end

function h=showvalues(obj,ax,x,y,val,space)
if ~obj.getSingleGuiParameter('showtext')
    h=[];
    return
end
if nargin<6
    space='';
end
for k=length(val):-1:1
    h(k)=text(ax,x(k),y(k),[space num2str(val(k),'%2.0f')],'FontSize',10,'Color','magenta','HitTest','off');
end
end

function selectrange(a,b,obj)
if strcmp(obj.axstep.Parent.Parent.SelectedTab.Title,'xy')
    rangex=obj.axxy.XLim;
    x=obj.coord.xr;
    indmin=find(x-obj.xyplotoffset<rangex(1),1,'last')+1;
    if isempty(indmin)
        indmin=1;
    end
    indmax=find(x-obj.xyplotoffset>rangex(2),1,'first')-1;
    if isempty(indmax)
        indmax=length(x);
    end
    t=obj.coord.timeplot;
    obj.range=[t(indmin) t(indmax)];
else
    obj.range=obj.axstep.XLim;
end
refit(a,b,obj,1)
 plotstatistics(obj)
end

% function [istepfit,svalfit]=splitmergefit(x,stepsize,p,steps)
% if contains(p.fitmode.selection,'mean')
%     mfun=@mean;
% else
%     mfun=@median;
% end
% 
% istep=[1 ; steps.properties.IndexStep];
% sval=[steps.properties.LevelBefore; steps.properties.LevelAfter(end)];
%     svalfit=sval;
%     istepfit=istep;
% for s=1:10
%     [istepfit, svalfit]=mergesplit(istepfit,svalfit,stepsize);
%     % recalculate sval based on x and istep2
%     [svalfit, istepfit]=fitstepind(x,istepfit,mfun);
%    
% end
% % 
% % if nargin>3 %y passed on
% %     svalfity=stepvalue(y,istepfit,mfun);
% % else 
% %     svalfity=[];
% % end
% end

% function [svalfit, istepfit]=fitstepind(x,istepfit,mfun)
%     for k=1:10
%         svalfit=stepvalue(x,istepfit,mfun);
%         istepfitold=istepfit;
%         istepfit=moveind(x,svalfit,istepfit);   
%         if all(istepfit==istepfitold)
%             k
%             break
%         end
%     end
% end

% function sval=stepvalue(x,istep,fun)
% istep=[istep ;length(x)];
% for k=length(istep)-1:-1:1
%     sval(k,1)=fun(x(istep(k)+1:istep(k+1)));
% end
% end

% function istep=moveind(x,sval,istep)
% istep=[istep; length(x)];
% % d=2;
% 
% % smin=0;splus=0;
% for k=2:length(istep)-1
%     dmin=0;
%     dplus=0;
%     smin=inf;splus=inf;
%     ih=istep(k);
% 
%     if abs(x(ih)-sval(k))<abs(x(ih)-sval(k-1)) %ih>1 && abs(x(ih-1)-sval(k))<abs(x(ih-1)-sval(k-1))
%         smin=abs(x(ih)-sval(k));
%         dmin=-1;
%     elseif ih>1 && abs(x(ih-1)-sval(k))<abs(x(ih-1)-sval(k-1)) %ih>1 && abs(x(ih-1)-sval(k))<abs(x(ih-1)-sval(k-1))
%         dmin=-2;
%         smin=abs(x(ih-1)-sval(k));
% 
% % %         smin=abs(x(ih-1)-sval(k))-abs(x(ih-1)-sval(k-1));
% %     elseif ih>2 && abs(mean(x(ih-2:ih-1)-sval(k)))<abs(mean(x(ih-2:ih-1)-sval(k-1)))
% %         dmin=-2;
% %         smin=abs(mean(x(ih-2:ih-1)-sval(k)))-abs(mean(x(ih-2:ih-1)-sval(k-1)));
%     end
%     if ih<length(x) && abs(x(ih+1)-sval(k-1))<abs(x(ih+1)-sval(k))% ih<length(x) && abs(x(ih)-sval(k))<abs(x(ih)-sval(k-1))
%         dplus=1;
%         splus=abs(x(ih+1)-sval(k-1));
%         
%     elseif ih<length(x)-1 && abs(x(ih+2)-sval(k-1))<abs(x(ih+2)-sval(k))
%         splus=abs(x(ih+2)-sval(k-1));
%         dplus=2;
% %         splus=abs(x(ih)-sval(k))-abs(x(ih)-sval(k-1));
% %     elseif ih<length(x)-1 && abs(mean(x(ih:ih+1)-sval(k)))<abs(mean(x(ih:ih+1)-sval(k-1)))
% %         dplus=2;
% %         splus=abs(mean(x(ih+1:ih+2)-sval(k)))-abs(mean(x(ih+1:ih+2)-sval(k-1)));
%     end
% %         ind=find(abs(x(max(1,ih-d):ih-1)-sval(k))<abs(x(max(1,ih-d):ih-1)-sval(k-1)));
% %         if ~isempty(ind)
% %             dmin=-d-1+min(ind);
% %         end
% %         lx=length(x);
% %         ind=find(abs(x(ih+1:min(lx,ih+d))-sval(k))>abs(x(ih+1:min(lx,ih+d))-sval(k-1)));
% %         if ~isempty(ind)
% %            dplus=max(ind);
% %         end
%      if smin<splus
%          istep(k)=max(1,max(istep(k)+dmin, istep(k-1)+1));
%      else 
%          istep(k)=max(1,min(istep(k)+dplus,istep(k+1)-1));
%      end
% end
% istep=istep(1:end-1);
% end

% function [istep, sval]=mergesplit(istep,sval,stepsize)
% stepv=diff(sval);
% step2=stepv(1:end-1)+stepv(2:end);
% sstep=find(abs(step2)<stepsize*1.4 & abs(stepv(1:end-1))<stepsize*.7);
% k=1;
% sstepc=sstep;
% while k<length(sstepc)
%     ind=find(sstepc==sstepc(k)+1);
%     if ~isempty(ind)
%         sstepc(ind)=[];
%         %k=k+1;
%     end
%     k=k+1;
% end
% istep(sstepc+2)=round((istep(sstepc+1)+istep(sstepc+2))/2);
% istep(sstepc+1)=[];
% sval(sstepc+1)=[];
% bigstep=find(abs(diff(sval))>stepsize*1.4)+1; %later pass on min / max step size
% 
% for k=1:length(bigstep)
%     [sval,istep]=insertstep(sval,istep,bigstep(k));
%     bigstep=bigstep+1; % as inserted, correct indices
% end
% end

% function [sval,istep]=insertstep(sval,istep,insertind)
%     istep=[istep(1:insertind) ;istep(insertind) ;istep(insertind+1:end)];
%     sval=[sval(1:insertind-1) ;(sval(insertind-1)+sval(insertind))/2 ;sval(insertind:end)];  
%     istep(insertind)=istep(insertind)-1;
%     istep(insertind+1)=istep(insertind+1)+1;
% end
% 
function manualcurate(obj)

f=figure(234);
ax=gca;
if ~isempty(obj.manualcuration)
    obj.manualcuration=ax.XLim;
end

hold(ax,'off')
if strcmp(obj.getSingleGuiParameter('filtertrackmode').selection,'smooth stepfind')
    xp=obj.coord.xr(obj.coord.indtime);
else
    xp=obj.coord.xfr(obj.coord.indtime);
end
plot(ax,obj.coord.timeplot(obj.coord.indtime),xp,'k','HitTest','off')
hold(ax,'on')

if isempty(obj.manualcuration)
    obj.manualcuration=ax.XLim;
end
ax.XLim=obj.manualcuration;

stairs(ax,obj.steps.steptime,obj.steps.stepvalue,'r','LineWidth',2,'HitTest','off')
stepv2=[obj.steps.stepvalue(1);(obj.steps.stepvalue(2:end)+obj.steps.stepvalue(1:end-1))/2];
showvalues(obj,ax,obj.steps.steptime(2:end),stepv2(2:end),diff(obj.steps.stepvalue),'    ');
roi=images.roi.Polyline(ax,'Position',horzcat(obj.steps.steptime,stepv2),'LineWidth',1,'Color','y');
addlistener(roi,'ROIMoved',@obj.updateroiposition);
addlistener(roi,'VertexAdded',@obj.addvertexroi);
addlistener(roi,'VertexDeleted',@obj.deletevertexroi);
plotsteps(obj)
end


function [sval,istep]=removestep(sval,istep,insertind)
if insertind+1<=length(istep)
    istep(insertind+1)=round((istep(insertind)+istep(insertind+1))/2);
end
    istep(insertind)=[];
    sval(insertind)=[];
end

function splitmerge(a,b,obj,what)
obj.manualcuration=[];
manualcurate(obj)
return

dcm_obj = datacursormode(obj.axstep.Parent.Parent.Parent);
info=dcm_obj.getCursorInfo;
if isempty(info)
    warndlg('Select step with the data tip tool first')
    return
end

steps=obj.steps;
[~,insertind]=min(abs(steps.steptime-info.Position(1)));
istep=obj.steps.indstep;
switch what
    case 1 %split
        [sval,istep]=insertstep(obj.steps.stepvalue,obj.steps.indstep,insertind);
    case 2
        [sval,istep]=removestep(obj.steps.stepvalue,obj.steps.indstep,insertind);
    case 3 %left       
        istep(insertind)=istep(insertind)-1;
    case 4 %right
        istep(insertind)=istep(insertind)+1;
end
if what<=2
    [svalfit, istep]=fitstepind(obj.coord.xr(obj.coord.indtime),istep,str2func(obj.getSingleGuiParameter('fitmode').selection));
end
calculatestepparameters(obj,istep);
plotsteps(obj);
end

function steps=calculatestepparameters(obj,stepindex)
if  ~isempty(obj.range)
    indtime=obj.coord.timeplot>=obj.range(1) & obj.coord.timeplot<=obj.range(2);
else
    indtime=true(size(obj.coord.xr)); 
end
obj.coord.indtime=indtime;
mfun=str2func(obj.getSingleGuiParameter('fitmode').selection);

if strcmp(obj.getSingleGuiParameter('filtertrackmode').selection,'smooth stepfind')
    x=obj.coord.xfr(indtime);
    y=obj.coord.yfr(indtime);
else
    x=obj.coord.xr(indtime);
    y=obj.coord.yr(indtime);
end



tv=obj.coord.timeplot(indtime);  
[stepv,nval]=stepvalue(x,stepindex,mfun);
steps.indstep=stepindex;
fhi=find(obj.index);
fh=fhi(indtime);
steps.indstepglobal=fh(stepindex(2:end));
steps.allindices=fh;
steps.steptime=tv(stepindex);
steps.stepvalue=stepv;
steps.stepsize=diff(stepv);
steps.possteps.x=stepv;
steps.possteps.y=stepvalue(y,stepindex,mfun);
if ~isempty(obj.coord.z)
    z=obj.coord.z(indtime);
    steps.possteps.z=stepvalue(z,stepindex,mfun);
    steps.std.z=stepvalue(z,stepindex,@std);
    steps.stddetrend.z=stepvalue(diff(z),stepindex,@std)/sqrt(2);
end
steps.possteps.time=tv(stepindex);
steps.dwelltime=diff(steps.steptime);
if isfield(obj.site.evaluation,obj.name)
    out=obj.site.evaluation.(obj.name);
end

%calculate statistics:
steps.std.x=stepvalue(x,stepindex,@std);
steps.std.y=stepvalue(y,stepindex,@std);

steps.stddetrend.x=stepvalue(diff(x),stepindex,@std)/sqrt(2);
steps.stddetrend.y=stepvalue(diff(y),stepindex,@std)/sqrt(2);

% if ~isempty(obj.coord.xfr)
    xf=obj.coord.xfr(indtime);
    yf=obj.coord.yfr(indtime);

    steps.std.xf=stepvalue(xf,stepindex,@std);
    steps.std.yf=stepvalue(yf,stepindex,@std);
    
    steps.stddetrend.xf=stepvalue(diff(xf),stepindex,@std)/sqrt(2);
    steps.stddetrend.yf=stepvalue(diff(yf),stepindex,@std)/sqrt(2);
% end

%find outliers: bad time resolution
mint=2;
window=2;
badtime=steptimeresolution(diff(tv),stepindex,mint,window);

bi=find(badtime);
bi2=bi+badtime(bi);
biall=union(bi,bi2);
biall(biall<1)=[];
biall(biall>length(badtime))=[];
badsteps=false(size(badtime));
badsteps(biall)=true;
steps.badsteps=badsteps;

steps.numlocsstep=nval;
out.steps=steps;
obj.steps=steps;
out.statall=calculatestatistics(obj,obj.index);
out.stattrack=calculatestatistics(obj,obj.index,obj.coord.indtime);
out.tracklength=x(end)-x(1);
out.range=obj.range;
obj.site.evaluation.(obj.name)=out;
end

function out=steptimeresolution(t,stepindex,mint,window)
out=zeros(length(stepindex)-1,1);
for k=1:length(stepindex)-1
    if any(t(stepindex(k)+1-1:min(length(stepindex),stepindex(k)+1+window-1))>mint)
        out(k,1)=-1;
    end
    if any(t(max(1,stepindex(k+1)-window):stepindex(k+1))>mint)
        out(k,1)=+1;
    end
end
end

function makemovie(a,b,obj)
plotsimple=obj.getSingleGuiParameter('simplemovie');
indt=obj.coord.indtime;
time=obj.coord.timeplot(indt);
frametime=obj.getSingleGuiParameter('frametime');
if isempty(frametime)
    frametime=mode(diff(time));
end

%XXX filter
x=obj.coord.xr(indt);
y=obj.coord.yr(indt);

% % XXXXXX 
% nmax=500;
% nmin=10;
% x=x(nmin:nmax);y=y(nmin:nmax);time=time(nmin:nmax);
if obj.getSingleGuiParameter('filtertrackmode').Value>1  
    fmode=obj.getSingleGuiParameter('filtermode').selection;
    windowsize=obj.getSingleGuiParameter('filterwindow'); 
    xf=runningWindowAnalysis(time,x,time,windowsize,fmode); 
    yf=runningWindowAnalysis(time,y,time,windowsize,fmode); 
    linew=1;
    filtertrackmode=true;

%     timen=min(time):fw:max(time);
%     xx=bindata(time,x,timen,fmode);
%     yy=bindata(time,y,timen,fmode);
%     x=xf;y=yf;
else
    xf=x;yf=y;
    linew=0.5;
    filtertrackmode=false;
end


ts=min(time):frametime:max(time);
f=figure(99);
f.Position(1)=1;f.Position(3)=1280;
ax=gca;

delete(ax.Children)

axis(ax,'equal');
axis(ax,'ij');

xlim(ax,[min(x)-10 max(x)+10])
ylim(ax,[min(y)-10 max(y)+10])
hold(ax,'on')
plot(ax,[min(x)-5 min(x)+10-5], [max(y) max(y)],'k','LineWidth',3)
ax.XTick=[];
ax.YTick=[];
for k=1:length(ts)
    indh=time<=ts(k);
    xh=xf(indh);
    yh=yf(indh);
    th=time(indh);
    tpassed=ts(k)-ts(1);
    ht=text(ax,double(min(xf)),double(min(yf)),[num2str(tpassed,'%3.0f') ' ms'],'FontSize',15);
        
    indc=obj.steps.steptime<ts(k);
    cx=obj.steps.possteps.x(indc);
    cy=obj.steps.possteps.y(indc);
    
    if filtertrackmode
         xr=x(indh);
         yr=y(indh);
         hr=plot(ax,xr,yr,'Color',[1 1 1]*0.7,'LineWidth',.5);
         hold(ax,'on')
    end
    hd=plot(ax,xh(end),yh(end),'ro','MarkerFaceColor','r','MarkerSize',15);
    hold(ax,'on')
   
    if ~plotsimple
        hb=plot(ax,xh,yh,'bo','MarkerSize',5,'MarkerFaceColor','b'); 
        hc=plot(ax,cx,cy,'m+','MarkerSize',15,'LineWidth',6);
    end
    hl=plot(ax,xh,yh,'k','LineWidth',linew);
    drawnow
    Fr(k)=getframe(ax);
    delete(hd)
    delete(hl)
    delete(ht)
    if filtertrackmode
    delete(hr)
    end
    if ~plotsimple
        delete(hb)
        delete(hc)
    end
end
smlfile=obj.getPar('lastSMLFile');
if ~isempty(smlfile)
    pfad=fileparts(smlfile);
else
    pfad=fileparts(obj.locData.files.file(1).name);
end

[file,pfad]=uiputfile([pfad filesep '*.mp4']);
if file
    mysavemovie(Fr,[pfad  file],'FrameRate',30)
end 
end

function resetview(a,b,obj)
obj.steps_x_range=[];
obj.xy_range=[];
plotsteps(obj)
end

function pard=guidef(obj)
pard.linkt.object=struct('String','Link','Style','text');
pard.linkt.position=[1,3];
pard.link.object=struct('String',{{'group','id'}},'Style','popupmenu','Value',2);
pard.link.position=[1,3.5];
pard.link.Width=1.5;

pard.overshoott.object=struct('String','Coarsness','Style','text');
pard.overshoott.position=[2,1];
pard.overshoot.object=struct('String','.8','Style','edit');
pard.overshoot.position=[2,2];
pard.overshoot.Width=0.5;

pard.fitmodet.object=struct('String','fit using','Style','text');
pard.fitmodet.position=[2,2.5];
pard.fitmode.object=struct('String',{{'mean','median'}},'Style','popupmenu');
pard.fitmode.position=[2,3.5];
pard.fitmode.Width=1.5;

pard.currentrange.object=struct('String','Current Range','Style','pushbutton','Callback',{{@selectrange,obj}});
pard.currentrange.position=[1,1];
pard.currentrange.Width=1.5;

pard.filterlocs.object=struct('String','filter (Renderer)','Style','checkbox');
pard.filterlocs.position=[4,1];
pard.filterlocs.Width=2;

pard.refit.object=struct('String','Refit','Style','pushbutton','Callback',{{@refit,obj,1}});
pard.refit.position=[4,2.7];
pard.refit.Width=0.8;
pard.refitalways.object=struct('String','always refit','Style','checkbox');
pard.refitalways.position=[4,3.5];
pard.refitalways.Width=1.5;


pard.refine.object=struct('String','Refine','Style','pushbutton','Callback',{{@refit,obj,2}});
pard.refine.position=[5,1];
pard.refine.Width=1;

pard.showtext.object=struct('String','values','Style','checkbox');
pard.showtext.position=[5,3];
pard.showtext.Width=1;

pard.msdanalysis.object=struct('String','MSD','Style','checkbox');
pard.msdanalysis.position=[5,4];
pard.msdanalysis.Width=1;

%auto-fit
p(1).value=0; p(1).on={}; p(1).off={'splitmerget','splitmergestep'};
p(2).value=1; p(2).on=p(1).off; p(2).off={};
pard.splitmerge.object=struct('Value',0,'String','Split/merge','Style','checkbox','Callback',{{@obj.switchvisible,p}});
pard.splitmerge.position=[3,1];
pard.splitmerge.Width=1.5;
pard.splitmerget.object=struct('String','step','Style','text','Visible','off');
pard.splitmerget.position=[3,2.5];
pard.splitmergestep.object=struct('String','','Style','edit','Visible','off');
pard.splitmergestep.position=[3,3];
pard.splitmergestep.Width=0.5;

pard.split.object=struct('String','Manual','Style','pushbutton','Callback',{{@splitmerge,obj,1}});
pard.split.position=[3,3.5];
pard.split.Width=1.5;





% pard.merge.object=struct('String','Merge','Style','pushbutton','Callback',{{@splitmerge,obj,2}});
% pard.merge.position=[4,2.5];
% pard.merge.Width=1.5;
% 
% pard.left.object=struct('String','<-','Style','pushbutton','Callback',{{@splitmerge,obj,3}});
% pard.left.position=[4,4];
% pard.left.Width=0.5;
% 
% pard.right.object=struct('String','->','Style','pushbutton','Callback',{{@splitmerge,obj,4}});
% pard.right.position=[4,4.5];
% pard.right.Width=0.5;


pard.makemovie.object=struct('String','Make Movie','Style','pushbutton','Callback',{{@makemovie,obj}});
pard.makemovie.position=[7,1];
pard.makemovie.Width=1.5;

pard.frametimet.object=struct('String','frametime','Style','text');
pard.frametimet.position=[7,2.5];
pard.frametime.object=struct('String','1','Style','edit');
pard.frametime.position=[7,3.5];
pard.frametime.Width=0.5;

pard.simplemovie.object=struct('String','simple','Style','checkbox');
pard.simplemovie.position=[7,4];


p(1).value=1; p(1).on={}; p(1).off={'filterwindowt','filterwindow','filtermode'};
p(2).value=2; p(2).on=p(1).off; p(2).off={};
p(3)=p(2); p(3).value=3;
% pard.filtertrackmode.object=struct('String','Filter','Style','checkbox','Callback',{{@obj.switchvisible,p}});
% pard.filtertrackmode.position=[6,1];

pard.filtertrackmode.object=struct('String',{{'no smooth','smooth plot','smooth stepfind'}},'Style','popupmenu','Callback',{{@obj.switchvisible,p}},'Value',2);
pard.filtertrackmode.position=[6,1];
pard.filtertrackmode.Width=1.7;

pard.filterwindowt.object=struct('String','dt ms:','Style','text');
pard.filterwindowt.position=[6,2.7];
pard.filterwindowt.Width=0.7;
pard.filterwindow.object=struct('String','1','Style','edit');
pard.filterwindow.position=[6,3.3];
pard.filterwindow.Width=0.4;
pard.filtermode.object=struct('String',{{'mean','median'}},'Style','popupmenu');
pard.filtermode.position=[6,3.7];
pard.filtermode.Width=1.3;

pard.resetview.object=struct('String','Reset view','Style','pushbutton','Callback',{{@resetview,obj}});
pard.resetview.position=[8,1];
pard.resetview.Width=1.3;

pard.onlyvld.object=struct('String','only vld','Style','checkbox','Value',1);
pard.onlyvld.position=[8,3];
pard.onlyvld.Width=2;

pard.fromdc.object=struct('String','display only second color, ch:','Style','checkbox','Value',0);
pard.fromdc.position=[9,1];
pard.fromdc.Width=3;

pard.dcch.object=struct('String','1','Style','edit');
pard.dcch.position=[9,4];
pard.dcch.Width=1;

% pard.dxt.Width=3;
pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi','se_sitepixelsize'};
pard.plugininfo.type='ROI_Evaluate';
pard.plugininfo.description='';
end
