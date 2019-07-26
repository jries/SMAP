classdef fibrilStraightener<interfaces.SEEvaluationProcessor
    % This plug-in is dependent of the BALM_fibril_growth.
    % Green line is the original boundary
    % White line is the refined boundary

    properties
        boundary
    end
    methods
        function obj=fibrilStraightener(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function out=run(obj, inp)
            out=runFibrilStraightener(obj, inp);
        end
     
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end

end



function pard=guidef(obj)
    % pard.dxt.Width=3;
    pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi'};
    pard.plugininfo.type='ROI_Evaluate';

end