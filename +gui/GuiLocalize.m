classdef GuiLocalize<interfaces.GuiModuleInterface&interfaces.LocDataInterface
    properties
        mainworkflow
        customworkflows
        batchfile
        loadpar=true;
    end
    methods
        function obj=GuiLocalize(varargin)
            obj@interfaces.GuiModuleInterface(varargin{:})
        end
        function makeGui(obj)
            h.loctab=uitabgroup(obj.handle);
            obj.adjusttabgroup(h.loctab);

            if ispc
                l2=obj.guiPar.FieldHeight;
                h2=obj.guiPar.FieldHeight*.95;
                dhg=0;
                dhgx=3;
                 dh=56;
                 dhx=22;
            else
                l2=obj.guiPar.FieldHeight*1.2;
                h2=obj.guiPar.FieldHeight*1.2;
                 dh=62;
                 dhg=3;
                 dhx=0;
                 dhgx=0;
            end
                       
            h.loctab.Units='pixel';
            h.loctab.Position(2)=h.loctab.Position(2)+dh;
            h.loctab.Position(4)=h.loctab.Position(4)-dh;
            h.loctab.Units='normalized';
            h.previewbutton=uicontrol(obj.handle,'Style','pushbutton','String','Preview','Position',[10 2, 70 h2],...
                'FontSize',obj.guiPar.fontsize,'Callback',{@preview_callback,obj},'FontWeight','bold');
            h.previewbutton.TooltipString=sprintf('Test current fitter settings on a single frame. \n Opens window which shows found localizations.\n Use before running fit on all frames.\n Use to open image to select ROIs for fitting.'); 
            h.previewframe=uicontrol(obj.handle,'Style','edit','String','1','Position',[160, 2, 60 h2],...
                'FontSize',obj.guiPar.fontsize,'Callback',{@previewframe_callback,obj});
           h.previewframeslider=uicontrol(obj.handle,'Style','slider','Position',[80, 2, 80, 22],...
                'FontSize',obj.guiPar.fontsize,'Callback',{@previewframeslider_callback,obj});
            
            h.localizebutton=uicontrol(obj.handle,'Style','pushbutton','String','Localize','Position',[250 2, 100, h2],...
                'FontSize',obj.guiPar.fontsize,'Callback',{@localize_callback,obj},'FontWeight','bold');
            h.localizebutton.TooltipString=sprintf('Start localization.');
            h.batchbutton=uicontrol(obj.handle,'Style','pushbutton','String','Batch','Position',[370+dhx, 2, 64, h2],...
                'FontSize',obj.guiPar.fontsize,'Callback',{@batch_callback,obj});
            h.previewframe.TooltipString=sprintf('Frame for preview. Select directly or with slider.');
            h.previewframeslider.TooltipString=h.previewframe.TooltipString;
            h.batchbutton.TooltipString=sprintf('Save batch file with all fitting parameters.');            
            
            h.batchprocessor=uicontrol(obj.handle,'Style','pushbutton','String','Processor','Position',[430+dhx, 2, 90 h2],...
                'FontSize',obj.guiPar.fontsize,'Callback',{@batchprocessor_callback,obj});
            h.batchprocessor.TooltipString=sprintf('Open the batch processor GUI to fit many files automatically with pre-defined settings.');
            h.wfinfo=uicontrol(obj.handle,'Style','pushbutton','String','Info','Position',[450+dhx, l2,50, h2],...
                'FontSize',obj.guiPar.fontsize,'Callback',{@wfinfo_callback,obj});
            h.wfinfo.TooltipString=sprintf('Show information about current workflow.');  
            
            h.wfcontext=uicontextmenu(getParentFigure(obj.handle));
            uimenu(h.wfcontext,'label','Info on workflow','Callback',{@wfinfo_callback,obj});
            uimenu(h.wfcontext,'label','Change workflow','Callback',{@wfload_callback,obj});
            uimenu(h.wfcontext,'label','Load raw WF','Callback',{@wfload_callback,obj});
            uimenu(h.wfcontext,'label','Load WF','Callback',{@wfload_callback,obj});
             uimenu(h.wfcontext,'label','Edit workflow','Callback',{@wfedit_callback,obj});
            uimenu(h.wfcontext,'label','Save workflow settings','Callback',{@savepar_callback,obj});
            uimenu(h.wfcontext,'label','Save workflow as...','Callback',{@wsave_callback,obj});
            h.wfname=uicontrol(obj.handle,'Style','text','String','x','Position',[10+10, l2, 350 h2],...
                'FontSize',obj.guiPar.fontsize,'HorizontalAlignment','left');    
            h.wfname.UIContextMenu=h.wfcontext;
            makemenuindicator(h.wfname,'lo');
            h.wfname.TooltipString='Name of current workflow. Right-click for context menu: info, change workflow, save current workflow settings';
             h.wfload=uicontrol(obj.handle,'Style','pushbutton','String','Change','Position',[370+dhx, l2, 80 h2],...
                'FontSize',obj.guiPar.fontsize,'Callback',{@wfload_callback,obj});
            %h.savepar=uicontrol(obj.handle,'Style','pushbutton','String','Save settings','Position',[290+dhx, l2, 80 h2],...
            %    'FontSize',obj.guiPar.fontsize,'Callback',{@savepar_callback,obj});
            
            h.wfsimple=uicontrol(obj.handle,'Style','togglebutton','String','-','Position',[500+dhx, l2+dhg,17+dhgx, h2-2*dhg],...
                'FontSize',obj.guiPar.fontsize,'Callback',{@wfsimplegui_callback,obj});
            h.wfsimple.TooltipString=sprintf('Show or hide advanced controls.');  
            
            outputfig=figure(207);
            outputfig.Visible='off';
            obj.setPar('loc_outputfig',outputfig)

            tabsizeh=obj.guiPar.tabsize2;
            tabsizeh(4)=tabsizeh(4)-dh;            
           
            settingsfile=obj.getGlobalSetting('mainLocalizeWFFile');
            par=readstruct(settingsfile);
            if ~isfield(par,'tab')
                par.tab=struct('hframe',struct('name','Input Image'),'hfilter',struct('name','Peak Finder'),'hfit',struct('name','Fitter'),'hloc',struct('name','Localizations'));
            end
            tabtags=fieldnames(par.tab);
            if isempty(tabtags)
                tabtags={'empty'};
                par.empty.name='empty';
            end
            if isfield(par,'workflowinfo')
                wfinfo=par.workflowinfo;
            else
                wfinfo='';
            end
            
            %create tabs
            for k=1:length(tabtags)
                h.(tabtags{k})=uitab(h.loctab,'Title',par.tab.(tabtags{k}).name);
                h.([tabtags{k} 'panel'])=uipanel(h.(tabtags{k}),'Unit','pixels','Position',tabsizeh); 
                replacestruct{k}={tabtags{k},h.([tabtags{k} 'panel'])};
            end
%             firstpanel=h.([tabtags{1} 'panel']);
             obj.guihandles=h;
            if ~exist(settingsfile,'file') %to make old settings work. Remove later.
                settingsfile=strrep(settingsfile,'settings',['settings' filesep 'workflows']);
            end
            par=readstruct(settingsfile,replacestruct);
            par=myrmfield(par,{'workflowinfo','tab'});
            if isempty(par)
                warndlg('cannot find settings file for fit workflow. Please set in menu SMAP/Preferences')
            else
                wffile=findsettingsfile(par.all.file,obj);
                [~,wfname]=fileparts(wffile);
                h.wfname.String=['Workflow: ' wfname];
           
                mainworkflow=interfaces.Workflow([],obj.P);
                mainworkflow.attachLocData(obj.locData);
                mainworkflow.setGuiAppearence(obj.guiPar)
                mainworkflow.setGuiAppearence(par.all)
                mainworkflow.processorgui=false;

    %             mainworkflow.pluginpath=wffile;
                mainworkflow.makeGui;
                mainworkflow.load(wffile,par,obj.loadpar);
                if ~isempty(wfinfo)
                    mainworkflow.description=wfinfo;
                end
                obj.mainworkflow=mainworkflow; 
                obj.children.mainworkflow=mainworkflow;
            end
            f=getParentFigure(obj.handle);
            c=uicontextmenu(f);
            h.loctab.UIContextMenu=c;
            m1 = uimenu(c,'Label','remove','Callback',{@menu_callback,obj});
            m3 = uimenu(c,'Label','add workflow','Callback',{@menu_callback,obj});
            
            previewframe_callback(0,0,obj)
            obj.addSynchronization('loc_fileinfo',[],[],@obj.update_slider)
            obj.addSynchronization('globalGuiState',[],[],@obj.update_guistate)

        end
        
        function update_slider(obj,a,b)        
            fi=obj.getPar('loc_fileinfo');
            if isempty(fi)
                return
            end
            nf=fi.numberOfFrames;
            if isempty(nf)||isinf(nf)||isnan(nf)
                nf=obj.getPar('loc_previewframe')+1;
%                   nf=2;
            end
            obj.guihandles.previewframeslider.Min=1;
            obj.guihandles.previewframeslider.Max=max(nf,1);
            pvf=max(1,min(nf,obj.getPar('loc_previewframe')));
            obj.guihandles.previewframeslider.Value=min(max(1,pvf),nf);
            obj.setPar('loc_previewframe',pvf)
            obj.guihandles.previewframe.String=num2str(pvf);
%             obj.setPar('loc_previewframe',round(pf));
            obj.guihandles.previewframeslider.SliderStep=[ceil(nf/50) ceil(nf/50)*10]/nf;
        end
        
        function setglobalguistate(obj,a,b)
            setglobalguistate@interfaces.GuiModuleInterface(obj,a,b);
            obj.guihandles.wfsimple.Value=obj.simplegui;
        end
        
        function update_guistate(obj)
            hs=obj.guihandles.wfsimple;
            switch obj.getPar('globalGuiState')
                case 's'
                    hs.Value=1;
                    hs.String='v';
                otherwise
                    hs.Value=0;
                    hs.String='-';
            end
                    
        end
    end
end
function wfinfo_callback(~,~,obj)
obj.mainworkflow.graph;
obj.mainworkflow.showinfo(true);
% if ~isempty(obj.mainworkflow.info.description)
% msgbox(obj.mainworkflow.info.description,obj.mainworkflow.info.name)
% end
end

function wfload_callback(a,b,obj)
if isprop(a,'Text') && contains(a.Text,'Load raw')
    obj.loadpar=false;
else
    obj.loadpar=true;
end

if (isprop(a,'Text') && contains(a.Text,'Change')) || (isprop(a,'String') && contains(a.String,'Change'))
    restorepar=true;
    oldsettings=obj.mainworkflow.getGuiParameters(true);
else
    restorepar=false;
end

obj.status('load workflow *.txt file');
drawnow
settingsfile=obj.getGlobalSetting('mainLocalizeWFFile');
[f,p]=uigetfile(settingsfile,'Select workflow *.txt file');
if f
    settingsfilen=[p f];
    loadwf(obj,settingsfilen)
    if restorepar
        obj.mainworkflow.setGuiParameters(oldsettings,true,false);
    end
    obj.status('new workflow loaded');
else
    obj.status('no workflow *.txt selected');
end

end

function wfedit_callback(~,~,obj)
 module=interfaces.Workflow;
module.processorgui=false;
module.handle=figure;
module.attachPar(obj.P);
module.attachLocData(obj.locData);

% p.Vrim=3;
% p.Xrim=4;
% p.FieldHeight=obj.guiPar.FieldHeight-4;
% module.setGuiAppearence(p)
module.makeGui;
wffile=obj.mainworkflow.pluginpath;
 module.load(wffile,[],false);
% obj.status('edit workflow *.txt file');
% drawnow
% settingsfile=obj.getGlobalSetting('mainLocalizeWFFile');
% [f,p]=uigetfile(settingsfile,'Select workflow *.txt file');
% if f
%     settingsfilen=[p filesep f];
%     loadwf(obj,settingsfilen)
%     obj.status('new workflow loaded');
% else
%     obj.status('no workflow *.txt selected');
% end

end


function loadwf(obj,settingsfilen)
    obj.setGlobalSetting('mainLocalizeWFFile',makerelativetopwr(settingsfilen));
    obj.mainworkflow.clear;
    delete(obj.handle.Children)
    obj.makeGui;
end
function menu_callback(callobj,~,obj)
switch callobj.Label
    case 'remove'
        selected=obj.guihandles.loctab.SelectedTab;
        title=selected.Title;
        if strcmp(title(1:2),'WF')
            number=str2double(title(3:end));
            delete(obj.customworkflows{number});
            obj.customworkflows(number)=[];
            delete(selected);
        else
            disp('only custom workflows are deletable')
        end

    case 'add workflow'
        nwf=length(obj.customworkflows);
        name=['WF' num2str(nwf+1)];
        obj.guihandles.(['tab_' name])=uitab(obj.guihandles.loctab,'Title',name);
        tabsizeh=obj.guiPar.tabsize2;
        tabsizeh(4)=235;
        obj.guihandles.(['panel_' name])=uipanel(obj.guihandles.(['tab_' name]),'Unit','pixels','Position',tabsizeh); 
        module=interfaces.Workflow;
        module.processorgui=false;
        module.handle=obj.guihandles.(['panel_' name]);
        module.attachPar(obj.P);
        module.attachLocData(obj.locData);
        p.Vrim=3;
        p.Xrim=4;
        p.FieldHeight=obj.guiPar.FieldHeight-4;
        module.setGuiAppearence(p)
        module.makeGui;
        obj.children.(name)=module;
        if nwf==0
            obj.customworkflows={module};
        else
            obj.customworkflows{end+1}=module;
        end
        
        obj.guihandles.loctab.SelectedTab= obj.guihandles.(['tab_' name]);
end
end

function preview_callback(a,b,obj)
selected=(obj.guihandles.loctab.SelectedTab.Title);
if strcmp(selected(1:2),'WF')
    number=str2double(selected(3:end));
    startmodule=obj.customworkflows{number};
else
    startmodule=obj.mainworkflow;
end
obj.setPar('loc_preview',true)
% obj.globpar.parameters.preview=true;
% notify(obj.P,'loc_initialize')
% startmodule.initialize;
obj.status('preview...');drawnow
startmodule.run;
obj.status('preview done');
end

function localize_callback(a,b,obj)
selected=(obj.guihandles.loctab.SelectedTab.Title);
if strcmp(selected(1:2),'WF')
    number=str2double(selected(3:end));
    startmodule=obj.customworkflows{number};
else
    startmodule=obj.mainworkflow;
end
obj.setPar('loc_preview',false)
% obj.globpar.parameters.preview=false;
% notify(obj.P,'loc_initialize')
% startmodule.initialize;
startmodule.run;
if ~isempty(obj.locData.loc)
obj.locData.regroup;
obj.setPar('locFields',fieldnames(obj.locData.loc));
obj.status('fitting done');
maingui=obj.getPar('mainGui');
maingui.setmaintab(3);
end
end

function batch_callback(a,b,obj)

[f,p]=uiputfile([fileparts(obj.getPar('filelist_localize')) '_batch.mat']);
if ~f
    return;
end

if strcmp(obj.guihandles.loctab.SelectedTab.Title,'Workflow')
    wf=obj.customworkflow;
else
    wf=obj.mainworkflow;
end

wf.save([p f])
obj.batchfile=[p f];
end

function previewframe_callback(a,b,obj)
pf=str2double(obj.guihandles.previewframe.String);
obj.setPar('loc_previewframe',round(pf));
if pf>1&&pf<obj.guihandles.previewframeslider.Max
obj.guihandles.previewframeslider.Value=round(pf);
end
end

function previewframeslider_callback(a,b,obj)
obj.guihandles.previewframe.String=num2str(round(obj.guihandles.previewframeslider.Value));
obj.setPar('loc_previewframe',round(obj.guihandles.previewframeslider.Value));
end
% 
% function plugino=makeplugin(obj,pluginpath,handle,pgui)
% % p=obj.guiPar;
% % p=copyfields(p,pgui);
% plugino=plugin(pluginpath{:});
% plugino.attachPar(obj.P);
% 
% % plugin.addGuiParameters(p);
% plugino.attachLocData(obj.locData);
% 
% plugino.setGuiAppearence(pgui);
% plugino.handle=handle;
% plugino.makeGui;
% 
% pluginname=pluginpath{end};
% obj.children.(pluginname)=plugino;
% end


function wfsimplegui_callback(object,~,obj)
if object.Value
    object.String='v';
else
    object.String='-';
end
obj.mainworkflow.fieldvisibility('guistate',object.Value)

end
function batchprocessor_callback(a,b,obj)
batchprocessor=WorkflowModules.Batchprocessor([],obj.P);
if ~isempty(obj.batchfile)
    batchprocessor.mainbatchfile=obj.batchfile;
end
batchprocessor.attachLocData(obj.locData);
batchprocessor.makeGui;
end

function savepar_callback(a,b,obj)
wf=obj.mainworkflow;
wf.save;
end

function wsave_callback(a,b,obj)
settingsfile=obj.getGlobalSetting('mainLocalizeWFFile');
[f,p]=uiputfile(settingsfile,'Select workflow *.txt file');
if ~f
    return
end
newtext=[p f];
[~,newname]=fileparts(f);
newmat=[p newname '.mat'];
psettings=readstruct(settingsfile);
psettings.all.file=newmat;
writestruct(newtext,psettings);

wf=obj.mainworkflow;
% wfmat=wf.pluginpath;
wf.save(newmat)
loadwf(obj,newtext)
obj.status(['workflow saved as ' newmat]);
end