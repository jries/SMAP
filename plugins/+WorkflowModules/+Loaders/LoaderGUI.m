classdef LoaderGUI<interfaces.WorkflowModule
    % LoaderGUI
    % Intermediate GUI to select a loader'
    properties
        loaders
        currentloader
        loadernames={'TifLoader','TifLoaderParallel','GrabFijiStacks','SimulateCameraImages'};
    end
    methods
        function obj=LoaderGUI(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
            obj.inputChannels=1; %1: image. 2: background image
            obj.isstartmodule=true;
            obj.excludeFromSave={'loaderlist'};
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
            obj.loadloaders;
        end
        function prerun(obj,p)
            p=obj.getGuiParameters;
            p=copyfields(p,obj.currentloader.getAllParameters);
            obj.currentloader=obj.loaders{p.loaderlist.Value};
            obj.currentloader.outputModules=obj.outputModules;
            obj.currentloader.prerun(p);    
        end
        function output=run(obj,data,p)
            p=copyfields(p,obj.currentloader.getGuiParameters);
            output=[];
            obj.currentloader.run(data,p);
        end
        function  addFile(obj,file)
             obj.currentloader.addFile(file);
        end
        function setoutputfilename(obj)
            obj.currentloader.setoutputfilename;
        end
        
        function loadloaders(obj)
            t1=obj.plugininfo.description;
            loadernames=obj.loadernames;
            obj.guihandles.loaderlist.String=obj.loadernames;
            p=obj.guiPar;
            p.Vrim=00;
            p.Vpos=1;
            poslist=obj.guihandles.loaderlist.Position;
            hpos=obj.handle.Position;
            panelpos(1)=poslist(1);
            panelpos(2)=poslist(2)-poslist(4)*6.5;
            panelpos(3)=hpos(3)-poslist(1)*2;
            panelpos(4)=poslist(4)*6.5;
            ip={};
            for k=1:length(loadernames)
                loader=plugin('WorkflowModules','Loaders',loadernames{k});
                loader.attachPar(obj.P);
                loader.attachLocData(obj.locData);
                hp=uipanel(obj.handle,'Units','pixels','Position',panelpos,'Visible','off');
                hp.Units='normalized';
                loader.handle=hp;
                loader.setGuiAppearence(p);
                loader.simplegui=obj.simplegui;
                loader.makeGui;
                ip=horzcat(ip,loader.inputParameters);
                obj.children.(loadernames{k})=loader;
                obj.loaders{k}=loader;
                obj.guihandles.([loadernames{k} '_panel'])=hp;
                t1=[t1 13 '   ' 96+k '. ' loader.info.name loader.info.description];
            end
            obj.inputParameters=ip;
            obj.loaders{1}.handle.Visible='on';
            obj.currentloader=obj.loaders{1};
            obj.plugininfo.description=t1;
        end
        function fieldvisibility(obj,varargin)
            fieldvisibility@interfaces.GuiModuleInterface(obj,varargin{:});
            fn=fieldnames(obj.children);   
            for k=1:length(fn)
                    obj.children.(fn{k}).fieldvisibility(varargin{:});
            end
        end
    end
end

function loaderlist_callback(object,b,obj)
for k=1:length(obj.loaders)
    obj.loaders{k}.handle.Visible='off';
end
loader=obj.loaders{object.Value};
loader.handle.Visible='on';
obj.currentloader=loader;

end

function pard=guidef(obj)
pard.loaderlist.object=struct('Style','popupmenu','String',{{'a'}},'Callback',{{@loaderlist_callback,obj}});
pard.loaderlist.position=[1,1];
pard.loaderlist.Height=1;
pard.loaderlist.Width=2;

pard.loaderlist.TooltipString=sprintf('Select here which loader plugin is used.');

pard.outputParameters={'loc_fitOnBackground'};
pard.plugininfo.type='WorkflowModule';
t1='Intermediate GUI to select a loader.';
pard.plugininfo.description=t1;
end