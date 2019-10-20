classdef Intensity2ManyChannels<interfaces.DialogProcessor
%     Assigns identity (color) to localizations based on intensities in
%     both camerachannels
    properties
        rois
        axis
        roihandles
        imparameters
        loadfile
    end
    methods
        function obj=Intensity2ManyChannels(varargin)        
            obj@interfaces.DialogProcessor(varargin{:}) ;  
            obj.history=true;
            obj.showresults=true;
            obj.propertiesToSave={'rois'};
        end
        function out=run(obj,p)
            out=[];
            if obj.processorgui==false || p.assignfield1.Value==1%run from WF
                p.assignfield1.selection='fit_n2';
                p.assignfield2.selection='fit_n1';
            end
            field1=p.assignfield1.selection;
            field2=p.assignfield2.selection;
            dx=obj.imparameters.dx;
            Nmax=double(ceil(obj.imparameters.max/dx));
            if p.usegrouped
                n1r=obj.locData.grouploc.(field1);
                n2r=obj.locData.grouploc.(field2);
            else
                n1r=obj.locData.loc.(field1);
                n2r=obj.locData.loc.(field2);
            end
            n1r(n1r<1)=1;
            n2r(n2r<1)=1;
            n1=log10(n1r)/dx;
            n2=log10(n2r)/dx;
            n1(n1<1)=1;
            n2(n2<1)=1;
            n1(n1>Nmax)=Nmax;
            n2(n2>Nmax)=Nmax;

            roi_callback(0,0,obj,0);
            channel=n1*0;
            for k=1:length(obj.rois)
                if ~isempty(obj.rois{k})
                    xpol=obj.rois{k}.Position(:,1)/dx;
                    ypol=obj.rois{k}.Position(:,2)/dx;
                    
                    imbw=poly2mask(xpol,ypol,Nmax,Nmax);
                    linind=sub2ind(size(imbw),round(n1),round(n2));
                    ischannel=imbw(linind);
                    channel(ischannel)=k;
                end
            end
            if p.usegrouped %convert back to single
                obj.locData.grouploc.channel=channel;
                channelu= obj.locData.grouped2ungrouped(1:length(channel),channel);
                obj.locData.loc.channel=channelu;
            else
                obj.locData.loc.channel=channel;
                obj.locData.regroup;
%                 obj.locData.grouploc.channel=channel(obj.locData.grouploc.groupindex);
            end
            
%             
            obj.setPar('locFields',fieldnames(obj.locData.loc))
%             
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function initGui(obj)
            obj.addSynchronization('locFields',[],[],@obj.updateLocFields)
            sdir=obj.getPar('SettingsDirectory');
            fout=[sdir filesep 'temp' filesep 'AssignColorRois.mat'];
            if exist(fout,'file')
                obj.loadfile=fout;
            end
%             obj.addSynchronization('filelist_short',obj.guihandles.dataselect,'String')
            obj.updateLocFields;
        end

        function updateLocFields(obj)
            excludefields={'frame','x','xnm','y','ynm','z','znm','locprecxnm','locprecynm',...
                'locprecznm','channel','PSFxnm','PSFynm','peakfindxnm','peakfindynm','filenumber',...
                'groupindex','numberInGroup','locprecnm'};
            if ~isempty(obj.locData.loc)
                fnpresent=fieldnames(obj.locData.loc);
                showfields=setdiff(fnpresent,excludefields);
                obj.guihandles.assignfield1.String=showfields;
                obj.guihandles.assignfield2.String=showfields;
                obj.guihandles.assignfield1.Value=min(obj.guihandles.assignfield1.Value,length(obj.guihandles.assignfield1.String));
                obj.guihandles.assignfield2.Value=min(obj.guihandles.assignfield2.Value,length(obj.guihandles.assignfield2.String));
                if ~isempty(obj.loadfile)
                    load_callback(0,0,obj,obj.loadfile)
                end
            end
        end

    end
end

function save_callback(a,b,obj)
sdir=obj.getPar('SettingsDirectory');
fout=[sdir filesep 'temp' filesep 'AssignColorRois.mat'];
[f,p]=uiputfile(fout);
rois=obj.rois;
field1=obj.getSingleGuiParameter('assignfield1').selection;
field2=obj.getSingleGuiParameter('assignfield2').selection;
if f
    save([p f],'rois','field1','field2')
end
end

function load_callback(a,b,obj,fout)
if nargin<4
sdir=obj.getPar('SettingsDirectory');
fout=[sdir filesep 'temp' filesep 'AssignColorRois.mat'];
[f,p]=uigetfile(fout);

if f
    fout=[p f];
    obj.loadfile=fout;
else
    return;
end
end
l=load(fout);
% obj.loaded=l;
obj.rois=l.rois;
field1=obj.getSingleGuiParameter('assignfield1');
f1v=find(strcmp(field1.String,l.field1));
if ~isempty(f1v)
    field1.Value=f1v(1);
    obj.setGuiParameters(struct('assignfield1',field1));
end
field2=obj.getSingleGuiParameter('assignfield1');
f2v=find(strcmp(field2.String,l.field2));
if ~isempty(f2v)
    field2.Value=f2v(1);
    obj.setGuiParameters(struct('assignfield2',field2));
end

if nargin<4 %only when called from load.
roi_callback(0,0,obj,0);
end

end

function roi_callback(a,b,obj,channel)

%channel==0; only draw

[n1,n2]=getintensities(obj);
% f=figure(944);
% ax=gca;
ax=obj.initaxis('n1 vs n2')
obj.axis=ax;
drawhistogram(obj,ax,n1,n2,obj.getSingleGuiParameter('logscale'))
cmap=prism(100);

if channel==-1
    obj.rois={};
    return
end

for k=1:length(obj.rois)
    if ~isempty(obj.rois{k}) 
        hroi=images.roi.Polygon(ax,'Position',obj.rois{k}.Position,'Color',cmap(k+1,:));
%         hroi.Tag=num2str(k);
        addlistener(hroi,'MovingROI',@(src,evt) updatePosition(src,evt,obj,k));
        obj.roihandles{k}=hroi;
    end
end
if channel>0
if length(obj.rois)<channel || isempty(obj.rois{channel})
    hroi=images.roi.Polygon(ax,'Color',cmap(channel+1,:));
    addlistener(hroi,'MovingROI',@(src,evt) updatePosition(src,evt,obj,channel));
%     hroi.Tag=num2str(channel);
    draw(hroi);
    obj.rois{channel}.Position=hroi.Position;
    
    obj.roihandles{channel}=hroi;
end
end
end

function updatePosition(src,evt,obj,channel)
% channel=str2double(src.Tag);
obj.rois{channel}.Position=src.Position;
end

function drawhistogram(obj,ax,n1,n2,islog)
mx=(max(quantile(n1,.99995), quantile(n2,.99995)));
% m0=(min(quantile(n1,.005), quantile(n2,.005)));

dp=0.02;
m0=-dp;
obj.imparameters.min=m0;
obj.imparameters.max=mx;
obj.imparameters.dx=dp;

range=m0:dp:mx;

him=histcounts2(n1,n2,range,range);
if islog
    him=log(him);
else
%     him=(him);
    q=quantile(him(:),0.9995);
    him(him>q)=q;
end

imagesc(ax,range,range,him)

end
function [n1,n2]=getintensities(obj)
field1=obj.getSingleGuiParameter('assignfield1').selection;
field2=obj.getSingleGuiParameter('assignfield2').selection;

[ll,ind]=obj.locData.getloc({'xnm',field1,field2},'Position','roi','layer',find(obj.getPar('sr_layerson')),'removeFilter','channel','grouping',obj.getSingleGuiParameter('usegrouped'));
n1r=(ll.(field1));
n2r=(ll.(field2));
n1r(n1r<1)=1;
n2r(n2r<1)=1;
n1=log10(n1r);
n2=log10(n2r);
end

function pard=guidef(obj)

pard.ch1t.object=struct('String','Channel 1','Style','text');
pard.ch1t.position=[1,3];
pard.ch1roi.object=struct('String','ROI 1','Style','pushbutton','Callback',{{@roi_callback,obj,1}});
pard.ch1roi.position=[1,4];

pard.ch2t.object=struct('String','Channel 2','Style','text');
pard.ch2t.position=[2,3];
pard.ch2roi.object=struct('String','ROI 2','Style','pushbutton','Callback',{{@roi_callback,obj,2}});
pard.ch2roi.position=[2,4];

pard.chNt.object=struct('String','Channel','Style','text');
pard.chNt.position=[3,3];
pard.chNt.Width=0.55;
pard.chN.object=struct('String','3','Style','edit');
pard.chN.position=[3,3.55];
pard.chN.Width=0.45;
pard.chNroi.object=struct('String','ROI N','Style','pushbutton','Callback',{{@roi_callback,obj,3}});
pard.chNroi.position=[3,4];

pard.showrois.object=struct('String','Show ROIs','Style','pushbutton','Callback',{{@roi_callback,obj,0}});
pard.showrois.position=[4,4];

pard.deleterois.object=struct('String','Delete ROIs','Style','pushbutton','Callback',{{@roi_callback,obj,-1}});
pard.deleterois.position=[5,4];

pard.t1.object=struct('String','fields','Style','text');
pard.t1.position=[1,1];

pard.assignfield1.object=struct('Style','popupmenu','String','n1|n2');
pard.assignfield1.position=[2,1];
pard.assignfield2.object=struct('Style','popupmenu','String','n1|n2');
pard.assignfield2.position=[3,1];

pard.logscale.object=struct('Style','checkbox','String','log scale','Value',1);
pard.logscale.position=[4,1];

pard.usegrouped.object=struct('Style','checkbox','String','use grouped','Value',0);
pard.usegrouped.position=[4,2];


pard.loadbutton.object=struct('String','load','Style','pushbutton','Callback',{{@load_callback,obj}});
pard.loadbutton.position=[5,1];
pard.savebutton.object=struct('String','save','Style','pushbutton','Callback',{{@save_callback,obj}});
pard.savebutton.position=[5,2];

pard.assignfield1.object.TooltipString='choose which field to use for splitting';
pard.assignfield2.object.TooltipString=pard.assignfield1.object.TooltipString;


pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.description='Assigns identity (color) to localizations based on intensities in both camerachannels';
end