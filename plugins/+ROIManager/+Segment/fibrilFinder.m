classdef fibrilFinder<interfaces.DialogProcessor&interfaces.SEProcessor
    % This is a plugin in development. Public has no access to the 
    % run_ functions called in this plugin. For internal users, 
    % "fibrilKymograph" is required. 
    methods
        function obj=fibrilFinder(varargin) 
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={};
        end
        function out=run(obj,p)
            out=runFibrilFinder(obj,p);
         
        end
        
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function pard = guidef(obj)
    pard.t1.object = struct('Style','text','String','test');
    pard.t1.position = [1,1];
    pard.t1.Width = 1;
    
    pard.plugininfo.type='ROI_Segment';
end