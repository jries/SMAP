classdef WFsplitter<interfaces.WorkflowModule
%     This plugin cuts out regions of interest of a defined size around the
%     candidate positions and passes these on to the fitter.
    properties
        loc_ROIsize
        preview
        disppreview
    end
    methods
        function obj=WFsplitter(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
            obj.inputChannels=2; 
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
   
            obj.setInputChannels(2,'frame');
            %populate here the splitWFselection menu. 
            %change visibility of sub-WF

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
    end
end



function pard=guidef
pard.splitWFselection.object=struct('Style','popupmenu','String',{{'1','2'}});
pard.splitWFselection.position=[1,1];
pard.splitWFselection.Width=0.8;



pard.syncParameters={{'splitWFselection','splitWFselection',{'Value'}}};

pard.plugininfo.type='WorkflowModule'; 
pard.plugininfo.description='This plugin cuts out regions of interest of a defined size around the candidate positions and passes these on to the fitter';
end