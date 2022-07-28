classdef GuiMainSMAP<interfaces.GuiModuleInterface & interfaces.LocDataInterface
%  DESCRIPTION:   SMAP: Main GUI
%  COPYRIGHT:     Jonas Ries, 2020
%  LICENSE:       GPLv3
%  AUTHOR:        Jonas Ries, EMBL Heidelberg, ries@embl.de 27.03.2020
%                 www.rieslab.de, www.github.com/jries/SMAP

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
            elseif ismac
                set(0,'DefaultUIControlFontSize',12);
            elseif isunix
                set(0,'DefaultUIControlFontSize',9);               
            end
            SMAP_stopnow=false;
            
            obj.setPar('maindirectory',pwd);
            disp(['main directory: ' pwd]);
             %settings directory 
             if isdeployed
                 if ispc
                    [status, result] = system('set PATH');
                    executableFolder = char(regexpi(result, 'Path=(.*?);', 'tokens', 'once'));
                    
                    programroot=ctfroot;
                    ind=strfind(programroot,filesep);
                    homedir=programroot(1:ind(3)-1);
                    possibledirs={[executableFolder filesep 'settings'],...
                        [ctfroot filesep 'settings'],...
                        [homedir filesep 'Documents' filesep 'settings'],...
                        [homedir filesep 'Documents' filesep 'MATLAB' filesep 'settings'],...
                        [homedir filesep 'Documents' filesep 'MATLAB' filesep 'SMAP' filesep 'settings'],...
                        [homedir filesep 'Documents' filesep 'SMAP' filesep 'settings'],...
                        'C:\Program Files\SMAP\application\settings',...
                        [homedir(1) ':\Program Files\SMAP\application\settings']};
            
                     
                 else
                     programroot=ctfroot;
                     ind=strfind(programroot,'application/SMAP.app/');
                     possibledirs={[programroot(1:ind-1) 'settings'],[programroot(1:ind-1) 'application' filesep 'settings']};
                 end
             else
                 if ispc
                      possibledirs={'settings',...
                     [fileparts(pwd) filesep  'settings'],...
                     [fileparts(fileparts(pwd)) filesep  'settings'],...
                     };
                 else
                     possibledirs={'settings',...
                     [pwd filesep 'MATLAB' filesep 'settings'],...
                     [pwd filesep 'MATLAB' filesep 'SMAP' filesep 'settings'],...
                     [pwd filesep 'Documents' filesep 'settings'],...
                     [pwd filesep 'Documents' filesep 'MATLAB' filesep 'settings'],...
                     [pwd filesep 'Documents' filesep 'MATLAB' filesep 'SMAP' filesep 'settings'],...
                     [pwd filesep 'Documents' filesep 'SMAP' filesep 'settings'],...
                     [pwd filesep  'settings'],...
                 };
                 end
             end
             settingsdir='';
             for k=1:length(possibledirs)
                 if exist(possibledirs{k},'dir')
                     settingsdir=possibledirs{k};
                     break
                 end
             end
             
             if isempty(settingsdir)
                 hwd=warndlg(['please select the directory /settings/ with the SMAP settings or copy the settings directory in any of the default directories and restart: ' possibledirs],'select settings','modal');
                 waitfor(hwd);
                 d=uigetdir(pwd,'select settings directory');
                 if d
                     settingsdir=d;
                 end
             end

             if k>1 %not default
                obj.setPar('maindirectory',fileparts(settingsdir));
             end
             disp(['settings directory: ' settingsdir]);

            settingsdirrel=makerelativetopwr(settingsdir);
            parentdir=fileparts(settingsdirrel);
            obj.setPar('SettingsDirectory',settingsdirrel);
            if ~isempty(parentdir)
                pluginhelp=[parentdir filesep 'Documentation' filesep 'help' ];
            else
                pluginhelp=[ 'Documentation' filesep 'help' ];
            end
            disp(['plugin help directory: ' pluginhelp]);
            obj.setPar('PluginHelpDirectory',pluginhelp);
                
            initglobalsettings(obj);
            if ~isdeployed
                addpath('shared');
                addpath('LocMoFit');
                addpath(pwd);
                if ~exist([settingsdir filesep 'temp'],'dir')
                    mkdir([settingsdir filesep 'temp'])
                end
                addpath('fit3Dcspline')
            else
                disp(pwd)
            end
        
            %update documentation from external files
            
            %from OC
%              urlzip='https://oc.embl.de/index.php/s/g0O4jQ4JEtmEris/download';
%              outdirdoc=[settingsdir filesep 'temp' filesep 'Documentation.tar'];
%              savewebfile(outdirdoc,urlzip);
%              if isdeployed
%                  outdir=[settingsdir filesep 'temp'];
%              else
%                  outdir=pwd;
%              end
%              unzip(outdirdoc,outdir);
%              delete(outdirdoc);
        %from tier1
            worked=false;
%            try 
                mainaddress='https://www.embl.de/download/ries/Documentation/';
                docfiles={'SMAP_manual_NPC.pdf','Example_SMAP_Step_by_step.pdf','ProgrammingGuide.pdf','SMAP_UserGuide.pdf','Getting_Started.pdf'};
%                 if isdeployed
%                      outdir=[settingsdir filesep 'temp' filesep 'Documentation' filesep];
%                 else
%                      outdir=[pwd filesep 'Documentation' filesep 'pdf' filesep];
%                 end
                 outdir= [fileparts(pluginhelp) filesep 'pdf' filesep];
                 if ~exist(outdir,'dir')
                      mkdir(outdir)
                 end
                for k=1:length(docfiles)
                    worked=worked|savewebfile([outdir docfiles{k}] ,[mainaddress docfiles{k}]);
                end
%             catch err
%                 err
%                 s
%             end
            if ~worked
                disp(['could not download and save documentation pdfs. Help might not work. Make sure you have write access to settings. Move the settings directory to ' possibledirs]);
                warndlg(['could not download and save documentation pdfs. Help might not work.  Make sure you have write access to settings. Move the settings directory to ' possibledirs])
            end
            
            %update plugin file if new plugins are saved
            makeplugincallfile('plugins');
            makeGeometricModelList;
            %add java path to bioformats
                bfpath=obj.getGlobalSetting('bioformatspath');
                bffile=[bfpath filesep 'bioformats_package.jar'];
                
            if exist(bffile,'file') %&& ~isdeployed
                javaaddpath(bffile);
            else
                disp('bioformats package not found. Please select path to bioformats_package.jar in the Preferences.')
                disp('you can download the Matlab toolbox for bioformats at  https://www.openmicroscopy.org/bio-formats/downloads/')
            end
           
            
            handle=figure(199);
             delete(handle.Children)
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
            tabpos=[2 32 obj.guiPar.width-2 368];

  
            gfile=obj.getGlobalSetting('guiPluginConfigFile');
            if ~exist(gfile,'file')
                ind = strfind(gfile,'settings');
                gfile=[settingsdir gfile(ind+8:end)];
%                 gfile=strrep(gfile,'settings', settingsdir);
            end
                
            if exist(gfile,'file')
                guimodules=readstruct(gfile,[],true);
            else
%                 guimodules=pmenu;
                guimodules=struct('File',[],'Analyze',[],'Process',[],'ROIManager',[]);
            end
            guimodulespure=myrmfield(guimodules,{'GuiParameters','globalGuiState'});
            obj.setPar('guimodules',guimodulespure);
            
            h.maintab = uitabgroup(handle,'Units','pixels','Position',tabpos);
            if ispc
                posmen='tri';
                shiftmen=[-10 -5];            
            elseif ismac
                posmen='tli';
                shiftmen=[10 -10];
            elseif isunix
                posmen='tri';
                shiftmen=[-10 -5];
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
            h.status.Tooltip='Progress of the current analysis';
            if ispc
                hstatus=30;
            else
                hstatus=31;
            end
            h.status.Position(4)=hstatus;
            h.status.Position(2)=1;
            obj.addSynchronization('status',h.status,{'String'})             

            
            h.errorindicator=uicontrol(handle,'Style','togglebutton','Units','normalized',...
                'Position',[0.01,0.002,.03,.03],'String',' ','Callback',{@error_reset,obj});
            h.errorindicator.Tooltip=sprintf('If an error occured during an analysis, this button turns \n red and you can read the error by clicking on it.');
            obj.addSynchronization('errorindicator',[],[],{@error_callback,obj,0}) 
            
            %Plugins

            obj.status('init plugins')
           
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

            gsites.maindir='ROIManager';
            gsites.guiplugins=myrmfield(guimodules.ROIManager,'Evaluate');
            gsites.makeGui();
            obj.children.guiSites=gsites;
            
            obj.guihandles=h;
            
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
            
            obj.status('all initialized')
            drawnow  
            set(handle, 'HandleVisibility', 'off');
            obj.setnormalizedpositionunits;
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

function error_callback(obj,b)
if strcmp(obj.getPar('errorindicator'),'clear')
    obj.guihandles.errorindicator.Value=0;
    resetstyle(obj)
else
    obj.guihandles.errorindicator.Value=1;
    obj.guihandles.errorindicator.BackgroundColor=[1 0 0];
    obj.guihandles.errorindicator.ForegroundColor=[1 0 0];
%     object.ForegroundColor=[0 1 1];
    obj.guihandles.errorindicator.FontWeight='bold';
    obj.guihandles.errorindicator.String='E';
end

end

function error_reset(a,b,obj)
if a.Value==0
warndlg(obj.getPar('errorindicator'));
resetstyle(obj)
end
end

function resetstyle(obj)
    obj.guihandles.errorindicator.BackgroundColor=[.94 .94 .94];
    obj.guihandles.errorindicator.ForegroundColor=[0 0 0];
    obj.guihandles.errorindicator.FontWeight='normal';
    obj.guihandles.errorindicator.String=' ';
end
% function saveplugin_callback(a,b,obj)
% 
% end