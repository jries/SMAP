classdef GuiPluginGroup< interfaces.GuiModuleInterface & interfaces.LocDataInterface
    properties
        maindir
%         plugindir
        resultspanel
        processors
        plugins
        guiplugins
        guipluginpath=['settings' filesep 'temp' filesep 'guimodules'];
    end

    methods
        function obj=GuiPluginGroup(varargin)
            obj@interfaces.GuiModuleInterface(varargin{:})
        end
        function makeGui(obj)
            obj.excludeFromSave(end+1)={'processorselect'};
            htab=obj.handle;
            set(htab,'Units','pixels')
                obj.guihandles.wrappanel=uipanel(obj.handle,'Unit','pixels','Position',obj.guiPar.tabsize2,'Visible','on');
                fontsize=obj.guiPar.fontsize;

                if ~isempty(obj.guiplugins)
                    %test for workflow : to plugin window
%                     if isfield(obj.guiplugins,'module')
%                        [~, name]=fileparts(obj.guiplugins.module);
%                         obj.makeplugin(name,{obj.guiplugins.module},0);
%                     else
                        
                        allplugins=fieldnames(obj.guiplugins);
                         allplugins=setdiff(allplugins,{'position','name','module'});
                         posm=[];
                        for k=1:length(allplugins)
                            if isfield(obj.guiplugins.(allplugins{k}),'position');
                                posm(k)=obj.guiplugins.(allplugins{k}).position;
                            else
                                posm(k)=inf;
                            end
    %                         if isfield(obj.guiplugins.(allplugins{k}),'name');
    %                             name{k}=obj.guiplugins.(allplugins{k}).name{1};
    %                         else
    %                             name{k}='';
    %                         end
                        end
                        [~,indpos]=sort(posm);
                        for k=1:length(allplugins)
                            obj.makeplugin(allplugins{indpos(k)});
                        end
%                     end
                    
                else
                    obj.plugins.pluginnames={'empty'};
                    obj.plugins.allclassnames={};
                end
                
                obj.initGui; 
                posparent=obj.guihandles.wrappanel.Position;
                obj.guihandles.processorselect=uicontrol(obj.guihandles.wrappanel,'Units','pixels',...
                    'Position',[5,posparent(4)-83,200,80],'Style','listbox','String',{'empty'},...
                    'FontSize',fontsize,'Callback',{@processorselect_callback,obj});
                
                if ~isempty(obj.guiplugins)&&~isempty(obj.plugins)
                    obj.setprocessorlist;
                    data.Value=1;
                    processorselect_callback(data,0,obj)
                end
                
            f=getParentFigure(obj.handle);
            c=uicontextmenu(f);
             obj.guihandles.processorselect.UIContextMenu=c;
             m4 = uimenu(c,'Label','move up','Callback',{@menu_callback,obj});
             m5 = uimenu(c,'Label','move down','Callback',{@menu_callback,obj});
             m1 = uimenu(c,'Label','add plugin','Callback',{@menu_callback,obj});
            m2 = uimenu(c,'Label','add workflow','Callback',{@menu_callback,obj});
            m6 = uimenu(c,'Label','rename','Callback',{@menu_callback,obj});
             m3 = uimenu(c,'Label','remove ','Callback',{@menu_callback,obj});
%              m4 = uimenu(c,'Label','reset: use all plugins','Callback',{@menu_callback,obj});
%              m8 = uimenu(c,'Label', 'load plugin structure','Callback',{@menu_callback,obj});
%              m9 = uimenu(c,'Label','save plugin structure','Callback',{@menu_callback,obj});
              m7 = uimenu(c,'Label','detach','Callback',{@menu_callback,obj});
            makemenuindicator(obj.guihandles.processorselect,'tri',[-15 -1]);
              
             
            
        end  
        function attachSE(obj,se)
            for k=1:length(obj.plugins.allclassnames)
                obj.processors.(obj.plugins.allclassnames{k}).attachSE(se);
            end
        end
        
        function thisplugin=makeplugin(obj,name,pluginpath,isprocessor)
%             if nargin<3
%                 screenname='';
%             end
%             obj.guiplugins.(name)
            if nargin<3
                pluginpath=obj.guiplugins.(name).module;
               
            end
            if nargin<4
                isprocessor=true;
            end
            if iscell(pluginpath)&&length(pluginpath)>=3
                thisplugin=[];
                try
                thisplugin=plugin(pluginpath{1:3});
                catch
                    
                    obj.guiplugins=myrmfield(obj.guiplugins,'name');
                end
                if isempty(thisplugin)||~isa(thisplugin,'interfaces.DialogProcessor')
                    if isstruct(thisplugin)
                        
                        obj.guiplugins=myrmfield(obj.guiplugins,'name');
                    else
                        disp([pluginpath{:} ' is not a subclass of interfaces.DialogProcessor'])
%                         pluginpath
%                         warning('selected plugin is not a subclass of interfaces.DialogProcessor')
                    end
                    return
                end
                isworkflow=false;
                processorgui=true;
            else %if iscell(pluginpath)&&length(pluginpath)==1
                thisplugin=interfaces.Workflow;
                thisplugin.pluginpath={pluginpath};
                
                isworkflow=true;
                processorgui=isprocessor;
            end
            
            
            panel=uipanel(obj.guihandles.wrappanel,'Unit','pixels','Position',obj.guiPar.tabsize2,'Visible','off');
            obj.guihandles.(name)=panel;
            
            thisplugin.attachLocData(obj.locData);
            thisplugin.attachPar(obj.P);
            

            thisplugin.handle=panel;
            parp.Vrim=85;
            thisplugin.setGuiAppearence(parp);
            thisplugin.processorgui=processorgui;
            thisplugin.makeGui;
            if isworkflow
                pluginpath
                thisplugin.load(pluginpath{1})
            end
            
            
            if isfield(obj.guiplugins.(name),'name')
                iname=obj.guiplugins.(name).name;
            else
            info=thisplugin.info;
            
            if ~isempty(info)
                iname=info.name;
            else
                iname=name;
            end
            end
            
            if isempty(obj.plugins)
                obj.plugins.allclassnames={name};       
                obj.plugins.pluginnames={iname};
            else
                numpl=length(obj.plugins.allclassnames);
                obj.plugins.allclassnames{numpl+1}=name;    
                obj.plugins.pluginnames{numpl+1}=iname;
            end
            
            obj.processors.(name)=thisplugin;
            obj.children.(name)=obj.processors.(name);
        end
        function setprocessorlist(obj)
%             if ~isempty(obj.plugins)
            obj.guihandles.processorselect.String=obj.plugins.pluginnames;
            if obj.guihandles.processorselect.Value<1||obj.guihandles.processorselect.Value>length(obj.plugins.pluginnames)
                obj.guihandles.processorselect.Value=1;
            end
%             end
        end
    end
end


function processorselect_callback(object,data,obj)
for k=1:length(obj.plugins.allclassnames)
    obj.processors.(obj.plugins.allclassnames{k}).setvisibility('off')
%     set(obj.guiPar.plugins.(obj.guiPar.allclassnames{k}).panel,'Visible','off')
end
% set(obj.guiPar.plugins.(obj.guiPar.allclassnames{object.Value}).panel,'Visible','on')
% obj.guiPar.plugins.(plugdir).(obj.guiPar.(plugdir).allclassnames{object.Value}).setvisibility('on')
obj.processors.(obj.plugins.allclassnames{object.Value}).setvisibility('on')
end

function menu_callback(callobj,b,obj)
guimodules=obj.getPar('guimodules');
% guimodules=readstruct('settings/temp/guimodules.txt',[],true);
switch callobj.Label
    case 'add plugin'
        plugins=obj.getPar('menu_plugins');
        types={'ProcessorPlugin','ROI_Analyze','SaverPlugin'};
        pg=browsefields(plugins,{},2,0,true,types);
        
        if isempty(pg)
            return
        end
        pluginname=field2cell(pg);
        name=pluginname{end};
        obj.guiplugins.(name).module=pluginname;
        obj.makeplugin(name);
        obj.setprocessorlist;
        lc=length(obj.guihandles.wrappanel.Children);
        obj.guihandles.wrappanel.Children(1:lc)=obj.guihandles.wrappanel.Children([2:lc 1]);
        guimodules.(obj.maindir{1}).(obj.maindir{2}).(name).module=pluginname;
        guimodules.(obj.maindir{1}).(obj.maindir{2}).(name).position=length(obj.plugins.pluginnames);
%         gm=load(obj.guipluginpath);
%         guimodules=gm.guimodules;
%         guimodules.(obj.maindir{1}).(obj.maindir{2}).(name)=pluginname;
        
%         guimodules.(pluginname{1}).(pluginname{2}).(pluginname{3})=pluginname;
%         save(obj.guipluginpath,'guimodules')
    case 'add workflow'
        [file,path]=uigetfile(['settings/workflows/*.mat']);
        if ~file
            return
        end
        [~,name]=fileparts(file);
        pluginname={[path filesep file]};
        obj.guiplugins.(name).module=pluginname;
        obj.makeplugin(name)
        obj.setprocessorlist;
        obj.guihandles.wrappanel.Children([1 2])=obj.guihandles.wrappanel.Children([2 1]);
        guimodules.(obj.maindir{1}).(obj.maindir{2}).(name).module=pluginname;
        guimodules.(obj.maindir{1}).(obj.maindir{2}).(name).position=length(obj.plugins.pluginnames);
        
    case 'remove '
        selection=obj.guihandles.processorselect.Value;
        name=obj.plugins.allclassnames{selection};
        delete(obj.processors.(name).handle);
%         ppath=obj.processors.(name).pluginpath;
         delete(obj.processors.(name));
        
        obj.plugins.allclassnames(selection)=[];
        obj.plugins.pluginnames(selection)=[];
        obj.guiplugins=rmfield(obj.guiplugins,name);
        obj.setprocessorlist;
        guimodules.(obj.maindir{1}).(obj.maindir{2})=rmfield(guimodules.(obj.maindir{1}).(obj.maindir{2}),name);
        obj.children=rmfield(obj.children,name);
        
    case {'move up','move down'}
        allnames=obj.plugins.allclassnames;
        selection=obj.guihandles.processorselect.Value;
%         name=allnames{selection};
        
        if strcmp(callobj.Label,'move up')
            newpos=max(1,selection-1);
        else
            newpos=min(length(allnames),selection+1);
        end
         obj.plugins.allclassnames([newpos selection])=obj.plugins.allclassnames([selection newpos]);
         obj.plugins.pluginnames([newpos selection])=obj.plugins.pluginnames([selection newpos]);
         obj.setprocessorlist;

    case 'rename'
        allnames=obj.plugins.pluginnames;
        selection=obj.guihandles.processorselect.Value;
        name=allnames{selection};
        newname=inputdlg('new name of module','rename module',1,{name});
        if isempty(newname)
            return
        end
        obj.plugins.pluginnames{selection}=newname{1};
        classname=obj.plugins.allclassnames{selection};
        guimodules.(obj.maindir{1}).(obj.maindir{2}).(classname).name=newname;
        obj.setprocessorlist;
%     case 'reset: use all plugins'
%         answ=questdlg('Delete all settings for pluginlist? Effect takes place after restarting the application');
%         
%         if strcmp(answ,'Yes')
%             obj.setGlobalSetting('guiPluginConfigFile','');
%         end
%             delete('settings/temp/guimodules.txt')
%             return
%         end
%     case 'load plugin structure'
%         fs='settings/*.txt';
%         [f,path]=uigetfile(fs,'load gui plugin structure file');
%         if f
%             guimodules=readstruct([path f],[],1);
%             msgbox('please restart SMAP')
%         end
%     case 'save plugin structure'
%         fs='settings/guiplugins.txt';
%         [f,path]=uiputfile(fs,'save gui plugin structure file');
%         if f
%              writestruct([path f],guimodules);
%         end
        
    case 'detach'
        f=figure('MenuBar','none','Toolbar','none');
        selection=obj.guihandles.processorselect.Value;
        module=obj.processors.(obj.plugins.allclassnames{selection});
        
        module.handle.Parent=f;

end
allclasses=obj.plugins.allclassnames;
for k=1:length(allclasses)
    guimodules.(obj.maindir{1}).(obj.maindir{2}).(allclasses{k}).position=k;
end
%  writestruct('settings/temp/guimodules.txt',guimodules);
 obj.setPar('guimodules',guimodules);
end


