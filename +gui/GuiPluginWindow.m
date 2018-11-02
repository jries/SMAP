classdef GuiPluginWindow< interfaces.GuiModuleInterface & interfaces.LocDataInterface
    properties
        guiplugins
         maindir
         tabgroup
%         plugindir
%         resultspanel
    end

    methods
        function obj=GuiPluginWindow(varargin)
            obj@interfaces.GuiModuleInterface(varargin{:})
        end
        
        function makeGui(obj)
            if isempty(obj.tabgroup)
                h.(obj.maindir)=uitabgroup(obj.handle);
            else
                h.(obj.maindir)=obj.tabgroup;
            end
            
%             htg=h.(obj.maindir);
            obj.adjusttabgroup(h.(obj.maindir));
%             if ispc
%             else
%                 htg.Units='pixel';
%                 htg.Position(1)=htg.Position(1)-8;
%                 htg.Position(3)=htg.Position(3)+16;
%                 htg.Position(2)=htg.Position(2)-12;
%                 htg.Position(4)=htg.Position(4)+16;
%                 htg.Units='normalized';
%             end
            obj.guihandles=h;
            guiplugins=obj.guiplugins;
            if ~isempty(guiplugins)
                fn=fieldnames(guiplugins);
                fn=setdiff(fn,{'position','name'});
                for k=1:length(fn)
                    if isfield(guiplugins.(fn{k}),'position')
                        pos(k)=guiplugins.(fn{k}).position;
                    else
                        pos(k)=inf;
                    end
                    if isfield(guiplugins.(fn{k}),'name')
                        name{k}=guiplugins.(fn{k}).name{1};
                    else
                        name{k}=fn{k};
                    end
                end
                [~,inds]=sort(pos);
                
                for k=1:length(fn)
                     obj.addplugingroup(fn{inds(k)},name{inds(k)})
                end
            end
                
           
            
            f=getParentFigure(obj.handle);
            c=uicontextmenu(f);
            h.(obj.maindir).UIContextMenu=c;
            m1 = uimenu(c,'Label','remove','Callback',{@menu_callback,obj});
            m2 = uimenu(c,'Label','add','Callback',{@menu_callback,obj});
            m3 = uimenu(c,'Label','add workflow','Callback',{@menu_callback,obj});
            m4 = uimenu(c,'Label','move left','Callback',{@menu_callback,obj});
            m5 = uimenu(c,'Label','move right','Callback',{@menu_callback,obj});
             m6 = uimenu(c,'Label','detach','Callback',{@menu_callback,obj});
            if ispc
                posmen='tri';
                shiftmen=[-10 -5];
            else
                posmen='tli';
                shiftmen=[10 -10];
            end
            makemenuindicator(h.(obj.maindir),posmen,shiftmen);
        end  
        function addplugingroup(obj,name,screenname)
            if nargin<3
                screenname=name;
            end
                    
                    ht=uitab(obj.guihandles.(obj.maindir),'Title',screenname,'Tag',name);
%                     ht.units='pixels';
%                     ht.Position(
%                     ht.units='normalized';
                    obj.guihandles.(['tab_' name])=ht;
%                     if isstruct(obj.guiplugins.(name))
                    if isstruct(obj.guiplugins.(name))&&~isfield(obj.guiplugins.(name),'module')    
                        pluging=gui.GuiPluginGroup(obj.guihandles.(['tab_' name]),obj.P);
                        pluging.guiplugins=obj.guiplugins.(name);
                        pluging.maindir={obj.maindir,name};
                        pluging.attachLocData(obj.locData);
%                         if ~isempty(obj.SE)
%                             pluging.attachSE(obj.SE)
%                         end
                        pluging.makeGui();
                    else
                        file=obj.guiplugins.(name);
                        pluging=makewf(obj,name,file);
                    end
                    obj.children.(name)=pluging;
            
        end
    end
end

function menu_callback(callobj,b,obj)
guimodules=obj.getPar('guimodules');
% guimodules=readstruct('settings/temp/guimodules.txt',[],true);
switch callobj.Label    
    case 'add'
        name=inputdlg('name of new tab');
        if isempty(name)
            return
        end
        name=name{1};
        obj.guiplugins.(name).position=length(fieldnames(obj.guiplugins))+1;
        obj.addplugingroup(name);
        guimodules.(obj.maindir).(name)=[];

    case 'remove'
        selected=obj.guihandles.(obj.maindir).SelectedTab;
        fieldr=selected.Title;
        if isfield(guimodules.(obj.maindir),fieldr) %only delete plugin groups
            delete(selected);
            guimodules.(obj.maindir)=rmfield(guimodules.(obj.maindir),fieldr); 
        end
    case 'add workflow'
        [file,path]=uigetfile(['settings/workflows/*.mat']);
        if file
            [~,name]=fileparts(file);
            name=name(1:min(end,16));
            obj.guiplugins.(name)=[path  file];
           addplugingroup(obj,name)
            guimodules.(obj.maindir).(name).module={[path  file]}; 
        else
            answ=inputdlg('Name of workflow','add empty workflow?');
            if ~isempty(answ)
                name=answ{1};
                obj.guiplugins.(name)=name;
                addplugingroup(obj,name);
                guimodules.(obj.maindir).(name).module={name}; 
            end
        end
    case {'move left','move right'}
        selected=obj.guihandles.(obj.maindir).SelectedTab;
        pos=find(selected.Parent.Children==selected);
        numtabs=length(selected.Parent.Children);
        if strcmp(callobj.Label,'move left')
            newpos=[1:pos-2,pos, pos-1,pos+1:numtabs];
        else
        
            newpos=[1:pos-1, pos+1,pos, pos+2:numtabs];
        end
        newpos(newpos>numtabs)=[];
        newpos(newpos<1)=[];
        selected.Parent.Children=selected.Parent.Children(newpos);
        
        for k=1:numtabs
            name=selected.Parent.Children(k).Title;
            guimodules.(obj.maindir).(name).position=k;
        end
    case 'detach'
        f=figure('MenuBar','none','Toolbar','none');
        hnew=uitabgroup(f);
        selected=obj.guihandles.(obj.maindir).SelectedTab;
        selected.Parent=hnew;
end
obj.setPar('guimodules',guimodules);
%  writestruct('settings/temp/guimodules.txt',guimodules);
end

function module=makewf(obj,name, file)
 
        
%         obj.guihandles.(['tab_' name])=uitab(obj.guihandles.(obj.maindir),'Title',name);
         module=interfaces.Workflow;
         module.processorgui=false;
         hpanel=uipanel(obj.guihandles.(['tab_' name]),'Position',[0 0 1 0.85]);
         module.handle=hpanel;%obj.guihandles.(['tab_' name]);
        module.attachPar(obj.P);
        module.attachLocData(obj.locData);
%         p.Vrim=40;
%         p.Xrim=4;
%         module.setGuiAppearence(p)
        module.makeGui;
        module.guihandles.showresults.Value=1;
        if exist(file,'file')
        module.load([file]);
        end
%         obj.children.(['tab_' name])=module;

end

