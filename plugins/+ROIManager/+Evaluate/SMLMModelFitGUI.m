classdef SMLMModelFitGUI<interfaces.SEEvaluationProcessor
    % This plug-in is dependent of the BALM_fibril_growth.
    % Green line is the original boundary
    % White line is the refined boundary

    properties
        boundary
    end
    methods
        function obj=SMLMModelFitGUI(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function out=run(obj, inp)
            out=runSMLMModelFitGUI(obj, inp);
        end
     
        function pard=guidef(obj)
            pard=runSMLMModelFitGUIdef(obj);
        end
    end

end