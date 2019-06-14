classdef GuiFile< interfaces.GuiModuleInterface & interfaces.LocDataInterface
    properties
        autosavetimer
        loaders
        savers
    end
    methods
        function obj=GuiFile(varargin)
            obj@interfaces.GuiModuleInterface(varargin{:})
            obj.outputParameters={'group_dx','group_dt'};
            obj.excludeFromSave={'filelist_long','loadmodule','savemodule'};
        end
        function initGui(obj)
            
            fw=obj.guiPar.FieldWidth;
            fh=obj.guiPar.FieldHeight;
            allloaders=pluginnames('File','Load');
            
            lp={};
            upperpos=obj.guihandles.load.Position;
            for k=1:length(allloaders)

                loaderhandle=uipanel(obj.handle,'Units','pixels','Position',[1,5.25*fh,2.9*fw,upperpos(2)-5.25*fh],'Visible','off');
                
                loaderp=plugin('File','Load',allloaders{k},loaderhandle,obj.P);
%                 infom.pluginpath={'File','Load',obj.loaders{k}};
                       
                loaderp.attachLocData(obj.locData);
                loaderp.makeGui;
                info=loaderp.info;     
                if ~strcmp(info.type,'LoaderPlugin')
                    delete(loaderhandle)
                else
                lp{k}=info.name;
                obj.guihandles.(['loaderpanel_' num2str(k)])=loaderhandle;
                obj.children.(['loader_' num2str(k)])=loaderp;
                loaders{k}=loaderp;
                end
            end
            obj.loaders=loaders;
            obj.loaders{1}.handle.Visible='on';
            obj.guihandles.loadmodule.String=lp;
            
            allsavers=pluginnames('File','Save');
            
            upperpos=obj.guihandles.save.Position;
%             obj.savers= plugintemp.plugins('File','Save',allsavers{1},[],obj.P); %to initialize
            for k=1:length(allsavers)
                saverhandle=uipanel(obj.handle,'Units','pixels','Position',[1,.6*fh,2.9*fw,upperpos(2)-.6*fh],'Visible','off');
               
                
                saver=plugin('File','Save',allsavers{k},saverhandle,obj.P);
%                 saver.pluginpath={'File','Save',allsavers{k}};
                saver.attachLocData(obj.locData);
                saver.makeGui;
                info=saver.info;
                
                if ~strcmp(info.type,'SaverPlugin')
                    delete(saverhandle)
                else
                    ls{k}=saver.info.name;
                     obj.guihandles.(['saverpanel_' num2str(k)])=saverhandle;
                    obj.children.(['saver_' num2str(k)])=saver;
                    savers{k}=saver;
                end
            end 
            obj.savers=savers;
            obj.savers{1}.handle.Visible='on';
            obj.guihandles.savemodule.String=ls;
            obj.guihandles.savemodule.Value=1;
            obj.addSynchronization('filelist_long',obj.guihandles.filelist_long,'String');
            obj.addSynchronization('autosavecheck',obj.guihandles.autosavecheck,'Value',{@autosavecheck_callback,0,0,obj});
            %make file menu
%             makefilemenu(obj);
            
            f=getParentFigure(obj.handle);
            c=uicontextmenu(f);
            obj.guihandles.filelist_long.UIContextMenu=c;
            obj.guihandles.filelist_long.String=[obj.getGlobalSetting('DataDirectory') filesep ];
            m1 = uimenu(c,'Label','info','Callback',{@menu_callback,obj});
            m1 = uimenu(c,'Label','remove','Callback',{@menu_callback,obj});
            m1 = uimenu(c,'Label','clear','Callback',{@menu_callback,obj});
            groupmode_callback(0,0,obj);
            
        end
     
        function loadbutton_callback(obj, handle,actiondata,isadd,pfad,f)
            %load data 
            p=obj.getAllParameters;
            fm=p.mainfile;   
            if isempty(fm)
                fm=p.filelist_long.selection;
            end
            path=fileparts(fm);          
            loader=obj.loaders{p.loadmodule.Value};
%             p=copyfields(p,loader.getGuiParameters);
%             loader=plugin('File','Load',obj.loaders{p.loadmodule.Value},[],obj.P);
%             loader.attachLocData(obj.locData);
%             loader.makeGui;
            try
                 ext=loader.info.extensions;
                 title=loader.info.dialogtitle;
            catch
                ext='*.*';
                title='format not specified';
                
            end
            if nargin<5
            [f,pfad]=uigetfile(ext,title,path,'MultiSelect','on');
            end
            if pfad %file selected
                if ~iscell(f)
                    f={f};
                end
                try
                 loader.clear([pfad f{1}],isadd)
                catch
                    obj.status('file type not recognized');
                    warning('file type not recognized');
%                     return
                end
%                 [~,~,ext]=fileparts(f{1});

%                 [mode, emptylocs]=getfilemode([pfad f{1}]);
%                 if isadd || ~emptylocs
%                     obj.locData.empty('filter');
%                 else
%                     obj.locData.empty;
%                 endgl
                
                for k=1:length(f)
                    [~,~,ext]=fileparts(f{k});
                    if ~isempty(strfind(f{k},'_sml'))||~isempty(strfind(f{k},'.csv'))
                        obj.setPar('mainfile',[pfad filesep f{k}]);
                    end
                    obj.status(['load: ' f{k}])
                    drawnow
%                     try
                    loadfiles(obj,loader,f{k},pfad)
%                     catch err
%                         err
%                         disp([f{k} ' could not be loaded'])
%                     end
                end

                obj.status('file loaded')
                if ~isempty(obj.locData.files.file)&&isfield(obj.locData.files.file(1),'transformation')
                    for k=length(obj.locData.files.file):-1:1
                        if ~isempty(obj.locData.files.file(k).transformation)
%                             obj.setPar('transformationfile','internal');
                            break
                        end
                    end
                end

                initGuiAfterLoad(obj)
                autosavecheck_callback(0,0,obj)
            end 
        end
        
       
        
        function remove_callback(obj,callobj,handle,actiondata)
            removefile=get(obj.guihandles.filelist_long,'Value');
            floc=obj.locData.loc.filenumber;
            removeind=floc==removefile;
            moveup=floc>removefile;
            
            obj.locData.loc.filenumber(moveup)= obj.locData.loc.filenumber(moveup)-1;
            obj.locData.removelocs(removeind);
            obj.locData.files.filenumberEnd=obj.locData.files.filenumberEnd-1;
            obj.locData.files.file(removefile)=[];
            obj.locData.SE.removeFile(removefile);
            
            fl=obj.getPar('filelist_long').String;
            fs=obj.getPar('filelist_short').String;
            fl(removefile)=[];
            fs(removefile)=[];
            obj.locData.regroup;
            obj.locData.filter;
            obj.setPar('filelist_long',fl,'String');
            obj.setPar('filelist_short',fs,'String');
            fsx=[{'layer','all'} fs];
            obj.setPar('filelist_short_ext',fsx,'String');
        end     
        
        function setGuiParameters(obj,p,varargin)
            setGuiParameters@interfaces.GuiModuleInterface(obj,p,varargin{:});
            if isfield(p,'filelist_long')
            obj.guihandles.filelist_long.String=p.filelist_long.String;   
            end
           
        end
        function group_callback(obj,a,b)
            obj.status('grouping');drawnow
            obj.locData.regroup;
            obj.status('grouping done');
%             group_callback(0, 0,obj);
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function delete(obj)
            delete(obj.autosavetimer)
        end
        
    end
end

function savemenu_callback(object,event,obj,what)
switch what
    case 'sml'
        savemodule=1;
    case 'tif'
        savemodule=2;
end
obj.guihandles.savemodule.Value=savemodule;
save_callback(0,0,obj)
end

function save_callback(object,event,obj)
p=obj.getAllParameters;
saver=obj.savers{p.savemodule.Value};
% psave=obj.getAllParameters(saver.inputParameters);
psave=saver.getAllParameters;
saver.save(psave);
end


function savemode_callback(data,b,obj)
for k=1:length(obj.savers)
    obj.savers{k}.handle.Visible='off';
end
obj.savers{data.Value}.handle.Visible='on';
end

function loadmode_callback(data,b,obj)
for k=1:length(obj.loaders)
    obj.loaders{k}.handle.Visible='off';
end
obj.loaders{data.Value}.handle.Visible='on';
end



function autosavecheck_callback(a,b,obj)

p=obj.getAllParameters;
t=obj.autosavetimer;
%creaste timer if empty or invalid
if isempty(t)||~isa(t,'timer')||~isvalid(t)
    t=timer;
    t.Period=p.autosavetime*60;
    t.StartDelay=t.Period;
    t.TimerFcn={@autosave_timer,obj};
    t.ExecutionMode='fixedRate';
    obj.autosavetimer=t;
end

if p.autosavecheck %deactivate
    disp('auto save enabled')
    if strcmpi(t.Running,'off')
    start(t);
    end
else
    if strcmpi(t.Running,'on')
        disp('auto save disabled')
    stop(t);
    end
end
end

function autosave_timer(a,b,obj)
p.mainGui=obj.getPar('mainGui');
p.saveroi=false;
if obj.guihandles.autosavecheck.Value
    savesml(obj.locData,'settings/temp/autosave_sml',p)
    time=datetime('now');
    disp(['autosave: ' num2str(time.Hour) ':' num2str(time.Minute)])
end
end

function autosavetime_callback(a,b,obj)
p=obj.getGuiParameters;
t=obj.autosavetimer;
if ~isempty(t)||isa(t,'timer')
    if strcmpi(t.Running,'on')
        stop(t);
    end
    obj.autosavetimer.Period=p.autosavetime*60;
    obj.autosavetimer.StartDelay=obj.autosavetimer.Period;
    if obj.guihandles.autosavecheck.Value
        start(t);
    end
end
end

function loadfiles(obj,loader,f,pfad)
mode=getfilemode([pfad f]);
if strcmp(mode,'tif')
    si=checkforsingleimages([pfad f]);             
    if si==1
        obj.setPar('filelist_localize',[pfad f]) %communication with localizer
        maing=obj.getPar('mainGui');
        maing.setmaintab(2);
        obj.locData.clear;
        return
    elseif si==0
%         return
    end
end


par=obj.getAllParameters(loader.inputParameters);
par=copyfields(par,loader.getGuiParameters);
loader.load(par,[pfad f]);
end

function menu_callback(menuobj,b,obj)
switch menuobj.Label
    case 'info'
        fin=obj.getPar('filelist_long').Value;
        info=obj.locData.files.file(fin).info;
        texta={};
        fn=fieldnames(info);
        for k=1:length(fn)
            ph=info.(fn{k});
            
            if iscell(ph)
                v='cell';
            elseif ischar(ph)
                v=ph;
            elseif isnumeric(ph)
                v=num2str(ph);
            end
                
            texta{end+1}=[fn{k} ': ' sprintf('\t') v];
        end
            
            listdlg('ListString',texta,'ListSize',[800,800]);
    case 'remove'
        obj.remove_callback;
    case 'clear'
        obj.locData.clear;
        fl={''};
        obj.setPar('filelist_long',fl,'String');
        obj.setPar('filelist_short',fl,'String');
        fsx=[{'layer','all'} fl];
        obj.setPar('filelist_short_ext',fsx,'String');

end
end

function makefilemenu(obj)
hsmap=obj.getPar('mainGuihandle');
hfile=uimenu(hsmap,'Label','File');
hload=uimenu(hfile,'Label','load','Callback',{@obj.loadbutton_callback,0});
hsavesml=uimenu(hfile,'Label','save SML','Callback',{@savemenu_callback,obj,'sml'});
hsavetif=uimenu(hfile,'Label','save Tif','Callback',{@savemenu_callback,obj,'tif'});

tobj=findobj(hsmap.Children,'Label','SMAP');
itarget=find(hsmap.Children==tobj);
ithis=find(hsmap.Children==hfile);
indold=1:length(hsmap.Children);
indnew=[ setdiff(indold,[ithis itarget])  ithis itarget];

% hdummy=hsmap.Children(itarget);
% hsmap.Children([itarget(1) ithis(1)])=hsmap.Children([ithis(1) itarget(1)]);
hsmap.Children=hsmap.Children(indnew);
% hsmap.Children(ithis)=hdummy;
end

function groupmode_callback(a,b,obj) 
% change visibility
% setPar(mode)
mode=obj.guihandles.group_mode.Value;
obj.setPar('group_mode',mode);

switch mode
    case 1
    obj.guihandles.group_dx.Visible='on';obj.guihandles.group_dxt.Visible='on';
    obj.guihandles.group_dxminmax.Visible='off';obj.guihandles.group_dxminmaxt.Visible='off';
    case 2
    obj.guihandles.group_dx.Visible='off';obj.guihandles.group_dxt.Visible='off';
    obj.guihandles.group_dxminmax.Visible='on';obj.guihandles.group_dxminmaxt.Visible='on';
end
end


function pard=guidef(obj)
pard.load.object=struct('Style','pushbutton','String','Load','FontWeight','bold','Callback',{{@obj.loadbutton_callback,0}});
pard.load.position=[4.75,2.5];
pard.load.Width=0.7;
pard.load.Height=1.5;

pard.add.object=struct('Style','pushbutton','String','Add','Callback',{{@obj.loadbutton_callback,1}});
pard.add.position=[4.75,3.2];
pard.add.Width=0.7;
pard.add.Height=1.5;

pard.loadmodule.object=struct('Style','popupmenu','String',{{'auto'}},'Callback',{{@loadmode_callback,obj}});
pard.loadmodule.position=[4.5,1];
 pard.loadmodule.Width=1.5;
 pard.loadmodule.object.TooltipString='select loader plugin';
 
% pard.updateGuiPar.object=struct('Style','checkbox','String','load Gui Parameters');
% pard.updateGuiPar.position=[5.5,1];
%  pard.updateGuiPar.Width=1.5;
% pard.updateGuiPar.TooltipString='Restore Gui parameters saved with localization data';

% 
% pard.remove.object=struct('Style','pushbutton','String','remove','Callback',{{@obj.remove_callback,1}});
% pard.remove.position=[4,4.5];
% pard.remove.Width=0.5;
% pard.add.Height=1.5;


pard.filelist_long.object=struct('Style','Listbox','String',{'x://'});
pard.filelist_long.position=[3,1];
pard.filelist_long.Width=4;
pard.filelist_long.Height=3;
if ispc
    psmen={'itr',[-20 -1]};
else
    psmen={'itr',[0 0]};
end
pard.filelist_long.uimenu=psmen;

% pard.uimenut.object=struct('Style','pushbutton','String','=');
% pard.uimenut.position=[1,1];
% pard.uimenut.Width=0.1;
% pard.uimenut.Height=.5;

pard.autosavecheck.object=struct('Style','checkbox','String','Auto save','Value',0);
pard.autosavecheck.position=[9,4];
pard.autosavecheck.Width=1;
pard.as_tdx.object=struct('Style','text','String','Interval min');
pard.as_tdx.position=[10,4];
pard.as_tdx.Width=0.65;
pard.autosavetime.object=struct('Style','edit','String','30','Callback',{{@autosavetime_callback,obj}});
pard.autosavetime.position=[10,4.65];
pard.autosavetime.Width=0.35;

pard.savemodule.object=struct('Style','popupmenu','String',{{'_sml','final image','raw images','_fitpos','settings'}},...
    'Callback',{{@savemode_callback,obj}});
pard.savemodule.position=[9,1.];
pard.savemodule.Width=1.5;

pard.save.object=struct('Style','pushbutton','String','Save','Callback',{{@save_callback,obj}});
pard.save.position=[9,2.5];
pard.save.Width=1;

pard.group_b.object=struct('Style','pushbutton','String','Group','Callback',{{@obj.group_callback}});
pard.group_b.position=[5,4];
pard.group_b.Width=1;

pard.group_tdx.object=struct('Style','text','String','dX');
pard.group_tdx.position=[7,4];
pard.group_tdx.Width=0.3;


pard.group_mode.object=struct('Style','popupmenu','String',{{'fix','locprec'}},'Callback',{{@groupmode_callback,obj}});
pard.group_mode.position=[7,4.2];
pard.group_mode.Width=0.8;

pard.group_dx.object=struct('Style','edit','String','35');
pard.group_dx.position=[8,4.3];
pard.group_dx.Width=0.35;
pard.group_dxt.object=struct('Style','text','String','nm');
pard.group_dxt.position=[8,4.65];
pard.group_dxt.Width=0.35;

pard.group_dxminmaxt.object=struct('Style','text','String','min max','Visible','off');
pard.group_dxminmaxt.position=[8,4];
pard.group_dxminmaxt.Width=0.5;

pard.group_dxminmax.object=struct('Style','edit','String','10 100','Visible','off');
pard.group_dxminmax.position=[8,4.5];
pard.group_dxminmax.Width=0.5;

pard.group_tdt.object=struct('Style','text','String','dT (frames)');
pard.group_tdt.position=[6,4];
pard.group_tdt.Width=0.65;
pard.group_dt.object=struct('Style','edit','String','1');
pard.group_dt.position=[6,4.65];
pard.group_dt.Width=0.35;


pard.syncParameters={{'filelist_long','filelist_long',{'String'}}};
pard.outputParameters= {'group_dx','group_dt','group_dxminmax','group_mode'};
pard.inputParameters={'mainfile'};


pard.load.object.TooltipString='load localization data or image. Load at least one localization data before adding a Tiff image.';
pard.add.object.TooltipString='add localization data or image';
% pard.remove.object.TooltipString='remove file';
pard.savemodule.object.TooltipString='select what to save';
pard.group_b.object.TooltipString='group localizations in consecutive frames';
pard.group_dx.object.TooltipString=sprintf('distance in nm which two locs can be apart \n and still grouped together. Empty: automatic');
pard.group_dt.object.TooltipString=sprintf('number of frames locs can be missing \n and still grouped together');

pard.group_mode.object.TooltipString=sprintf('Maximum allowed distance: 1. fix: set to fix value (below). \n 2. locprec: from localizaton precision. dx^2< 2.5(lp1^2+lp2^2) \n dx>min, dx<max');

pard.autosavecheck.object.TooltipString=sprintf('save localizations every XX minutes as settings/temp/autosave_sml.mat');
end