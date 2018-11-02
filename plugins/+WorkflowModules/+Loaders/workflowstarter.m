classdef workflowstarter<interfaces.WorkflowModule
    properties
    end
    methods
        function obj=workflowstarter(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
            obj.inputChannels=0;
            obj.isstartmodule=true;
        end

        function initGui(obj)

        end
        function prerun(obj,p)
           
        end
        function run(obj,data,p)
            dateof=interfaces.WorkflowData;
            dateof.frame=1;
            dateof.ID=1;
            dateof.eof=true;
            obj.output(dateof)
        end
        function pard=guidef(obj)
            pard.startb.object=struct('Style','pushbutton','String','Run','Callback',@obj.startcallback);
            pard.plugininfo.type='WorkflowModule'; 
            pard.startb.position=[1,1];
        end
        function startcallback(obj,a,b)
            obj.parent.run;
%              obj.initialize;
%             obj.run;
        end

    end
end