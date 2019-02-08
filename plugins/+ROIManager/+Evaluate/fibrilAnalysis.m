classdef fibrilAnalysis<interfaces.SEEvaluationProcessor
    % This is a plugin in development. Public has no access to the 
    % run_ functions called in this plugin. For internal users, 
    % "fibrilKymograph" is required.
    methods
        function obj=fibrilAnalysis(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        
        function out=run(obj,p)
            if p.filtering
                if isfield(obj.site.evaluation,'fibrilAnalysis')&&isfield(obj.site.evaluation.fibrilAnalysis,'setting')
                else
                    p.axisLb = 0;
                    p.axisUb = 0;
                end
            end
            out = runFibrilAnalysis(obj,p);
            
            if isfield(obj.site.evaluation,'fibrilAnalysis')&&isfield(obj.site.evaluation.fibrilAnalysis,'setting')
                obj.guihandles.axisLb.String = out.setting.axisLb;
                obj.guihandles.axisUb.String = out.setting.axisUb;
            else
                obj.guihandles.axisLb.String = 0;
                obj.guihandles.axisUb.String = out.straitened.arclength;
            end
                
            obj.guihandles.filtering.Value = 0;
            
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
    
    pard.t1.object = struct('Style','text','String','Filtering:');
    pard.t1.position = [2 1];
    pard.t1.Width = 1;
    
    pard.axisLb.object = struct('Style','edit','String','0');
    pard.axisLb.position = [2 2];
    pard.axisLb.Width = 0.5;
    
    pard.axisUb.object = struct('Style','edit','String','0');
    pard.axisUb.position = [2 2.5];
    pard.axisUb.Width = 0.5;
    
    pard.filtering.object = struct('Style','checkbox','String','','Value',0);
    pard.filtering.position = [2 3];
    pard.filtering.Width = 0.5;

    pard.redoProjection.object = struct('Style','checkbox','String','Re run','Value', 0);
    pard.redoProjection.position = [3 1];
    pard.redoProjection.Width = 1;

    pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi'};
    pard.plugininfo.type='ROI_Evaluate';
end