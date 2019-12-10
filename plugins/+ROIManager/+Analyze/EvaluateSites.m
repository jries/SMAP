classdef EvaluateSites<interfaces.DialogProcessor&interfaces.SEProcessor
%     Evaluates all sites. This is mainly used for the BatchAnalysis
%     plugin.
    methods
        function obj=EvaluateSites(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
        end
        
        function out=run(obj,p)  
            out=[];
           eval=obj.SE.processors.eval; 
           eval.redrawall;
           
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end




function pard=guidef

pard.t1.object=struct('String','run evaluation plugins for all ROIs','Style','text');
pard.t1.position=[1,1];
pard.t1.Width=1;


pard.plugininfo.type='ROI_Analyze';
pard.plugininfo.description='Evaluates all sites. This is mainly used for the BatchAnalysis plugin.';

end