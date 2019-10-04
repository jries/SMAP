classdef WFsplitter<interfaces.WorkflowModule
%     This plugin cuts out regions of interest of a defined size around the
%     candidate positions and passes these on to the fitter.
    properties
        modules
        guistates
    end
    methods
        function obj=WFsplitter(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
            obj.inputChannels=2; 
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
   
            obj.setInputChannels(2,'frame');
        end
        function prerun(obj,p)
         
        end
        function nooutput=run(obj,data,p)
           nooutput=[];
            emptydat=data{1};
            emptydat.data=[];
            switch p.splitWFselection.Value
                case 1
                    dato{1}=data{1};
                    dato{3}=data{2};
                    dato{2}=emptydat;
                    dato{4}=emptydat;
                    obj.output(dato{1},1);
                    obj.output(dato{3},3);
                case 2
                    dato{2}=data{1};
                    dato{4}=data{2};
                    dato{1}=emptydat;
                    dato{3}=emptydat;
                    obj.output(dato{2},2);
                    obj.output(dato{4},4);
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
                module=obj.outputModules(2*k).module;
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