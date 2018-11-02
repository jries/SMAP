classdef SEProcessor<interfaces.GuiModuleInterface & interfaces.LocDataInterface
    properties
        SE
    end
    methods
         function obj=SEProcessor(varargin)
            obj@interfaces.GuiModuleInterface(varargin{:})
        end
        function attachSE(obj,se)
            obj.SE=se;
            if isempty(obj.locData)
                obj.locData=obj.SE.locData;
            end
        end
        function attachLocData(obj,locData)
            attachLocData@interfaces.LocDataInterface(obj,locData);
            if isempty(obj.SE)
                obj.SE=locData.SE;
            end
        end
        function updateSingleParameter(obj, data,actionData,field)
            val=obj.getSingleGuiParameter(field);
            obj.SE.sePar.(data.Parent.Title).(field)=val;
        end
        function updateParameters(obj)
            fn=fieldnames(obj.guihandles);
            for k=1:length(fn)
                obj.updateSingleParameter(obj.guihandles.(fn{k}),0,fn{k})
            end           
        end

    end
end