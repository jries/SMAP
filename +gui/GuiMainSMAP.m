classdef GuiMainSMAP<interfaces.GuiModuleInterface & interfaces.LocDataInterface
    methods
        function obj=GuiMainSMAP(varargin)
            obj.attachPar(interfaces.ParameterData);
            obj.P.clear;
            obj.locData.attachPar(obj.P)
%             obj.locData.iscopy=false;

        end
        function makeGui(obj)
            global SMAP_stopnow
            if ispc
            set(0,'DefaultUIControlFontSize',9);
            else
                set(0,'DefaultUIControlFontSize',12);
            end
            SMAP_stopnow=false;
            addpath('shared');
            addpath(pwd);
            if ~exist(['settings' filesep 'temp'],'dir')
                mkdir(['settings' filesep 'temp'])
            end
            obj.setPar('maindirectory',pwd);
            
             %global settings
            initglobalsettings(obj);
            
            makeplugincallfile('plugins');
            
            %add java path
%             try
%                 bfpath='/Users/ries/Downloads/bfmatlab/';
                bfpath=obj.getGlobalSetting('bioformatspath');
            if exist(bfpath,'dir')
                addpath(bfpath)
                bfCheckJavaPath;
            else
%             catch err %no bioformats found
                disp('bioformats package not found. Please select path to bioformats_package.jar in the Preferences.')
                disp('you can download the Matalb toolbox for bioformats at  https://www.openmicroscopy.org/bio-formats/downloads/')
            end
           
            
            handle=figure(199);
%             clf(handle,'reset')
             delete(handle.Children)
%              

            obj.setPar('mainGuihandle',handle);
            obj.setPar('mainGui',obj);
            obj.setPar('synchronizeguistate',true);
            obj.handle=handle;
            set(handle,'MenuBar','none','Toolbar','none')
            [pmenu,hmenu]=makePluginMenu(obj,handle);
            obj.setPar('menu_plugins',pmenu);
            
            obj.guiPar.width=550;
            scrsz = get(groot,'ScreenSize');
            height=min(scrsz(4)-80,760);
            hpos=max(5,min(scrsz(4)-height-80,scrsz(4)/5));
            if ispc
                vpossmap=8;
            else
                vpossmap=3;
            end
            set(handle,'Units','Pixels', 'Name','SMAP')
            set(handle,'Position',[vpossmap hpos obj.guiPar.width height]);            
            set(handle,'ButtonDownFcn',{@figure_selected,obj},...
                'SizeChangedFcn',{@sizechanged_callback,obj},'NumberTitle','off')
drawnow
%             set(handle, 'Name','SMAP');
            tabpos=[2 32 obj.guiPar.width-2 368];

   
            
            gfile=obj.getGlobalSetting('guiPluginConfigFile');
            gfile=findsettingsfile(gfile);
            
            if exist(gfile,'file')
                guimodules=readstruct(gfile,[],true);

            else
                guimodules=pmenu;
%                 guimodules.globalGuiState='a';
            end
            guimodulespure=myrmfield(guimodules,{'GuiParameters','globalGuiState'});
            obj.setPar('guimodules',guimodulespure);
            
            
           
            
            
            h.maintab = uitabgroup(handle,'Units','pixels','Position',tabpos);
            if ispc
                posmen='tri';
                shiftmen=[-10 -5];
            else
                posmen='tli';
                shiftmen=[10 -10];
            end
            makemenuindicator(h.maintab,posmen,shiftmen);
            f=getParentFigure(obj.handle);
            ch=uicontextmenu(f);
            h.maintab.UIContextMenu=ch;
            m1 = uimenu(ch,'Label','detach','Callback',{@detach_callback,obj,h.maintab});
      

            h.tab_file = uitab(h.maintab,'Title','File');
            h.tab_loc = uitab(h.maintab,'Title','Localize');
            h.tab_render = uitab(h.maintab,'Title','Render');
            h.tab_process = uitab(h.maintab,'Title','Process','Tag','tab_process');
            h.tab_analyze = uitab(h.maintab,'Title','Analyze','Tag','tab_analyze');
            h.tab_siteexplorer=uitab(h.maintab,'Title','ROIs','Tag','tab_siteexplorer');
            set(h.maintab,'SelectedTab',h.tab_file)
            
            h.stopnow=uicontrol('Style','togglebutton','Units','normalized',...
                'Position',[0.9,0.002,.07,.03],'String','Stop','Callback',{@stopnow_callback,obj});
            h.stopnow.Units='pixels';
            h.stopnow.Position(4)=28;
            h.stopnow.TooltipString='Interrupt execution of certain commands (e.g. fitting). Not implemented for all plugins';
            h.status=uicontrol(handle,'Style','text','Units','normalized',...
                           'String','status','Position',[0 0 .8 0.035]);
            h.status.Units='pixels';
            if ispc
                hstatus=30;
            else
                hstatus=31;
            end
            h.status.Position(4)=hstatus;
            h.status.Position(2)=1;
            obj.addSynchronization('status',h.status,{'String'})             

            %Plugins

            obj.status('init plugins')
           
%             obj.locData.guiData.siteexplorer=gsites;
            
            
            %file
            h.filepanel=uipanel(h.tab_file,'Units','pixel','Position',obj.guiPar.tabsize1);
            obj.status('init GuiFile')
            gfile=gui.GuiFile(h.filepanel,obj.P);
            gfile.attachLocData(obj.locData);
            gfile.makeGui();
            obj.children.guiFile=gfile;
            
            %localize
            obj.status('init Localize')
            h.localizepanel=uipanel(h.tab_loc,'units','pixel','Position',obj.guiPar.tabsize1);
            gloc=gui.GuiLocalize(h.localizepanel,obj.P);
            gloc.attachLocData(obj.locData);
            gloc.makeGui();
            obj.children.guiLoc=gloc;
            
            %Render tab
            h.renderpanel=uipanel(h.tab_render,'units','pixel','Position',obj.guiPar.tabsize1);
            obj.status('init GuiRender')
            grec=gui.GuiRender(h.renderpanel,obj.P);
            grec.attachLocData(obj.locData);
            grec.makeGui();
            obj.children.guiRender=grec;
            
           
            
            %Process tab
            obj.status('init processor plugins')
            h.processpanel=uipanel(h.tab_process,'units','pixel','Position',obj.guiPar.tabsize1);
            gprocess=gui.GuiPluginWindow(h.processpanel,obj.P);
            gprocess.attachLocData(obj.locData);
            gprocess.guiplugins=guimodules.Process;
             gprocess.maindir='Process';
%             gprocess.plugindir={'drift','Register','Assign2C','Modify'};
            gprocess.makeGui();
            obj.children.guiProcess=gprocess;
            
%             %Analysis tab   
            obj.status('init analysis plugins')
            h.analysispanel=uipanel(h.tab_analyze,'units','pixel','Position',obj.guiPar.tabsize1);
            ganalysis=gui.GuiPluginWindow(h.analysispanel,obj.P);
            ganalysis.attachLocData(obj.locData);
            ganalysis.maindir='Analyze';
            ganalysis.guiplugins=guimodules.Analyze;
            ganalysis.makeGui();
            
            
            %statistics,xy vs time, blinking stats
            %line profile, image correlation stuff (line, image)
            %calibrate 3DA, side view, 3D volume
            % SPT, analyzePSF, movie

            obj.children.guiAnalysis=ganalysis;

            
             %Site explorer
            obj.status('init ROI manager')
            h.sitespanel=uipanel(h.tab_siteexplorer,'Units','pixel','Position',obj.guiPar.tabsize1);
            gsites=gui.SEMainGui(h.sitespanel,obj.P);
            gsites.attachLocData(obj.locData);
%             gsites.addGuiParameters(guipar);
            gsites.maindir='ROIManager';
            gsites.guiplugins=myrmfield(guimodules.ROIManager,'Evaluate');
            gsites.makeGui();
            obj.children.guiSites=gsites;
            
            obj.guihandles=h;
%             gui.setTooltips(obj);
            
            
            %undo
            undo=gui.Undo(obj.handle,obj.P);
            undo.attachLocData(obj.locData);
            undo.makeGui();
            obj.children.undo=undo;
            
             if isfield(guimodules,'GuiParameters')
                
                obj.setGuiParameters(guimodules.GuiParameters,true);
                
             end
            
            cb=hmenu.hsimplegui.Callback;
            if  isfield(guimodules,'globalGuiState')&&~isempty(guimodules.globalGuiState)&& strcmp(guimodules.globalGuiState,'s')
                feval(cb{1},struct('Checked','off'),0,cb{2})
                
            else
                feval(cb{1},struct('Checked','on'),0,cb{2})
            end
            
            
%             sizechanged_callback(obj.handle, 0, obj)
            obj.status('all initialized')
            drawnow  
            set(handle, 'HandleVisibility', 'off');
        end
        function psaved=saveParameters(obj)
            psaved=obj.getGuiParameters(true);
            if nargout==0
             save('parfile.mat','psaved');
            end
        end
        function loadParameters(obj)
            load parfile.mat
            obj.setGuiParameters(psaved,true);
        end
        function setmaintab(obj,number)
            tab=obj.guihandles.maintab.Children(number);
            obj.guihandles.maintab.SelectedTab=tab;
        end
        function delete(obj)
            disp('delete')
            t=timerfindall;
            if ~isempty(t)
                stop(t);
                delete(t);
            end
            delete(obj.getPar('sr_figurehandle'))
            delete(obj.handle)
        end
    end

end


function sizechanged_callback(object, event, obj)
uiwait(obj.handle,1)  
f=object.Position(3)/obj.guiPar.width;
if f~=1
obj.resize(f);
obj.guiPar.width=object.Position(3);
end
end

function figure_selected(object, event, obj)
figure(obj.getPar('sr_figurenumber'))
end

function stopnow_callback(object,b,obj)
global SMAP_stopnow
SMAP_stopnow=object.Value;
if object.Value
    object.BackgroundColor=[1 0 0];
%     object.ForegroundColor=[0 1 1];
    object.FontWeight='bold';
else
    object.BackgroundColor=[.94 .94 .94];
    object.ForegroundColor=[0 0 0];
     object.FontWeight='normal';
end

end

function detach_callback(a,b,obj,htab)
f=figure('MenuBar','none','Toolbar','none');
hnew=uitabgroup(f);
htab.SelectedTab.Parent=hnew;

end

% function saveplugin_callback(a,b,obj)
% 
% end