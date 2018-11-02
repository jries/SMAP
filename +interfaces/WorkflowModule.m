classdef WorkflowModule<interfaces.WorkflowInterface&interfaces.GuiModuleInterface&interfaces.LocDataInterface
    properties
        parent
        updatetime=1
        timervalue=tic;
    end
    methods
        function obj=WorkflowModule(varargin)
            obj@interfaces.GuiModuleInterface(varargin{:})
        end

        function initGui(obj)
            
        end
         function updateGui(obj)
         end

        function status(obj,str)
            if toc(obj.timervalue)>obj.updatetime
                status@interfaces.GuiModuleInterface(obj,str);
                drawnow limitrate;
                obj.timervalue=tic;
            end
        end
        function sethandle(obj,h)
            obj.handle=h;
        end

    end
end

