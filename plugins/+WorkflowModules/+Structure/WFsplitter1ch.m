classdef WFsplitter1ch<interfaces.WorkflowModule
%     This plugin cuts out regions of interest of a defined size around the
%     candidate positions and passes these on to the fitter.
    properties
        modules
        guistates
    end
    methods
        function obj=WFsplitter1ch(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
            obj.inputChannels=1; 
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
   
            obj.setInputChannels(1,'frame');
        end
        function prerun(obj,p)
         
        end
        function nooutput=run(obj,data,p)
           nooutput=[];
            switch p.splitWFselection.Value
                case 1
                    obj.output(data{1},1);
                case 2
                    obj.output(data{2},2);
            end
        end
        function modelchanged(obj,a,b)
            br=obj.getSingleGuiParameter('splitWFselection').Value;
            br2=2-br+1;
     
                for k=1:length(obj.modules{br})
                    guih=obj.modules{br}{k}.guihandles;                    
                    if isempty(guih)
                        continue
                    end

                    fn=fieldnames(guih);
                    for l=1:length(fn)
                        guih.(fn{l}).Visible='on';
                    end
                        obj.modules{br}{k}.switchvisibleall;
                end
                for k=1:length(obj.modules{br2})
                    guih=obj.modules{br2}{k}.guihandles;
                    if isempty(guih)
                        continue
                    end
                    fn=fieldnames(guih);
                    for l=1:length(fn)
                        guih.(fn{l}).Visible='off';
                    end
                end            
        end
        function updateGui(obj)
           %populate here the splitWFselection menu. 
            %change visibility of sub-WF
            for k=1:2
                module=obj.outputModules(k).module;
                ind=1;
                name{k}='';
                while ~contains(module.info.name,'WFcombiner')
                    obj.modules{k}{ind}=module;
                    name{k}=[name{k} ',' module.info.name];
                    ind=ind+1;
                    module=module.outputModules(1).module;
                end
                name{k}(1)='';
            end
           
            obj.setGuiParameters(struct('splitWFselection',struct('String',{name},'Value',1))  )  
            obj.modelchanged;
        end
    end
end



function pard=guidef(obj)
pard.splitWFselection.object=struct('Style','popupmenu','String',{{'1','2'}},'Callback',@obj.modelchanged);
pard.splitWFselection.position=[1,1];
pard.splitWFselection.Width=2;



pard.syncParameters={{'splitWFselection','splitWFselection',{'Value'}}};

pard.plugininfo.type='WorkflowModule'; 
pard.plugininfo.description='This plugin cuts out regions of interest of a defined size around the candidate positions and passes these on to the fitter';
end