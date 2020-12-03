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
            if isstruct(obj.plugininfo.description)
                fn=fieldnames(obj.plugininfo.description);
                t1=obj.plugininfo.description.(fn{1});
                fieldn=fn{1};
            else
                t1=obj.plugininfo.description;
                fieldn='none';
            end
            loadernames=obj.loadernames;
            obj.guihandles.loaderlist.String=obj.loadernames;
            p=obj.guiPar;
            p.Vrim=00;
            p.Vpos=1;
            poslist=obj.guihandles.loaderlist.Position;
            hpos=obj.handle.Position;
            panelpos(1)=poslist(1);
            panelpos(2)=poslist(2)-poslist(4)*5.5;
            panelpos(3)=hpos(3)-poslist(1)*2;
            panelpos(4)=poslist(4)*5.5;
            ip={};
            t1=[t1 '\n Specific loaders: '];
            for k=1:length(loadernames)
                loader=plugin('WorkflowModules','Loaders',loadernames{k});
                loader.attachPar(obj.P);
                loader.attachLocData(obj.locData);
                hp=uipanel(obj.handle,'Units','pixels','Position',panelpos,'Visible','off');
                
                loader.handle=hp;
                loader.setGuiAppearence(p);
                loader.simplegui=obj.simplegui;
                loader.makeGui;
                ip=horzcat(ip,loader.inputParameters);
                obj.children.(loadernames{k})=loader;
                obj.loaders{k}=loader;
%                 hp.Units='normalized';
                obj.guihandles.([loadernames{k} '_panel'])=hp;
                if isstruct(loader.info.description)
                    fn=fieldnames(loader.info.description);
                    th=loader.info.description.(fn{1});
                else
                    th=loader.info.description;
                end
                t1=[t1 '\n' 96+k '. ' loader.info.name '\n' th];
            end
            obj.inputParameters=ip;
            obj.loaders{1}.handle.Visible='on';
            obj.currentloader=obj.loaders{1};
%             obj.plugininfo
%             t1
%             obj.plugininfo.description
            obj.plugininfo.description.(fieldn)=t1;
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