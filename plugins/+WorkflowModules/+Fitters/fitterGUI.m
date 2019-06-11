classdef fitterGUI<interfaces.WorkflowModule
%     Intermediate GUI to select a fitting plugin.
    properties
        fitters
        currentfitter
    end
    methods
        function obj=fitterGUI(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
            obj.inputChannels=2; %1: image. 2: background image
            
            obj.outputParameters={'bg_dx','bg_dt','subtractbg'};
            obj.excludeFromSave={'fitterlist'};
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
            obj.guihandles.fitterlist.Callback={@fitterlist_callback,obj};
            obj.setInputChannels(obj.inputChannels,'frame');
            obj.loadfitters;
%            obj.guihandles.camparbutton.Callback={@camparbutton_callback,obj};
        end
        function prerun(obj,p)
            p=obj.getGuiParameters;
            obj.currentfitter=obj.fitters{p.fitterlist.Value};
            obj.currentfitter.outputModules=obj.outputModules;
            obj.currentfitter.prerun;  
        end
        function output=run(obj,data,p)
            output=[];
            obj.currentfitter.run(data);
        end

        function loadfitters(obj)
            t1=obj.plugininfo.description;
            fitnames={'MLE_GPU_Yiming','EMCCD_SE_MLE_GPU','RadialSymmetry2D','RadialSymmetry3D'};
            obj.guihandles.fitterlist.String=fitnames;
            p=obj.guiPar;
            p.Vrim=00;
            p.Vpos=1;
            poslist=obj.guihandles.fitterlist.Position;
            hpos=obj.handle.Position;
            panelpos(1)=poslist(1);
            panelpos(2)=0;
            panelpos(3)=hpos(3);
            panelpos(4)=poslist(2);
            for k=1:length(fitnames);
                fitter=plugin('WorkflowModules','Fitters',fitnames{k});
                fitter.attachPar(obj.P);
                hp=uipanel(obj.handle,'Units','pixels','Position',panelpos,'Visible','off');
                fitter.handle=hp;
                fitter.setGuiAppearence(p);
                fitter.simplegui=obj.simplegui;
                fitter.makeGui;
                obj.children.(fitnames{k})=fitter;
                obj.fitters{k}=fitter;
                obj.guihandles.([fitnames{k} '_panel'])=hp;
                t1=[t1 13 '   ' 96+k '. ' fitter.info.name fitter.info.description];
            end
            obj.fitters{1}.handle.Visible='on';
            obj.currentfitter=obj.fitters{1};
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

function fitterlist_callback(object,b,obj)
for k=1:length(obj.fitters)
    obj.fitters{k}.handle.Visible='off';
end
fitter=obj.fitters{object.Value};
fitter.handle.Visible='on';
obj.currentfitter=fitter;

end

function pard=guidef(obj)
pard.fitterlist.object=struct('Style','listbox','String','EMCCD: Single Emitter|sCMOS: Single Emitter|EMCCD: Multiple Emitter|sCMOS: Mulitple Emitter');
pard.fitterlist.position=[4,1];
pard.fitterlist.Height=3;
pard.fitterlist.Width=3;

pard.fitterlist.TooltipString=sprintf('Select here which Fitter plugin is used.');

pard.outputParameters={'loc_fitOnBackground'};
pard.plugininfo.type='WorkflowModule';
t1='Intermediate GUI to select a fitting plugin.';
pard.plugininfo.description=t1;
end