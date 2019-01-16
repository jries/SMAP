classdef SEProcessor<interfaces.GuiModuleInterface & interfaces.LocDataInterface
    properties
        SE
    end
    methods
         function obj=SEProcessor(varargin)
            obj@interfaces.GuiModuleInterface(varargin{:})
        end
        function attachSE(obj,se)
            obj.attachLocData(se.locData);
%             obj.SE=se;
%             if isempty(obj.locData)
%                 obj.locData.SE=obj.SE;
%             end
        end
        function attachLocData(obj,locData)
            
            attachLocData@interfaces.LocDataInterface(obj,locData);
            
            obj.locData.SE.attachPar(obj.locData.P);
%             if isempty(obj.SE)
%                 obj.SE=locData.SE;
%             end
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
        function out=get.SE(obj)
            out=obj.locData.SE;
        end
        function set.SE(obj,se)
            obj.locData.SE=se;
        end

    end
end