classdef usePeakfinder<interfaces.WorkflowModule
%    Passes on coordinates from peak finder
    properties
    end
    methods
        function obj=usePeakfinder(varargin)
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
        function outputdat=run(obj,data,p)
            outputdat=data{1};
           %do some conversions here
        end
    end
end


function pard=guidef(obj)

pard.text.object=struct('Style','text','String','use coordinates from peak finder (e.g. deepSMLM)');
pard.text.position=[1,1];
pard.text.Width=4;
pard.plugininfo.type='WorkflowModule'; 
pard.plugininfo.description='Passes on coordinates from peak finder';

end