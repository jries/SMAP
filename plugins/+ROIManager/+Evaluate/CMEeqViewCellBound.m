classdef CMEeqViewCellBound<interfaces.SEEvaluationProcessor
    properties

    end
    methods
        function obj=CMEeqViewCellBound(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function out=run(obj,p)
            out = runCMEeqViewCellBound(obj,p);
        end
        function pard=guidef(obj);
            pard=guidef;
        end
    end

end

function pard=guidef
pard.boun.object=struct('Style','checkbox','String','Show boundary');
pard.boun.position=[1,1];
pard.boun.Width=2;
pard.tan.object=struct('Style','checkbox','String','Show tangentLine');
pard.tan.position=[2,1];
pard.tan.Width=2;

% pard.dxt.Width=3;
pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi'};
pard.plugininfo.type='ROI_Evaluate';
end