classdef getKymogroph<interfaces.SEEvaluationProcessor
    % This is a plugin in development. Public has no access to the 
    % run_ functions called in this plugin. For internal users, 
    % "fibrilKymograph" is required.
    methods
        function obj=getKymogroph(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        
        function out=run(obj,p)
            out = runGetKymogroph(obj,p);
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function pard = guidef(obj)
    pard.showPlot.object = struct('Style','checkbox','String','Show fig', 'Value', 0);
    pard.showPlot.position = [1 1];
    pard.showPlot.Width = 1;

    pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi'};
    pard.plugininfo.type='ROI_Evaluate';
end