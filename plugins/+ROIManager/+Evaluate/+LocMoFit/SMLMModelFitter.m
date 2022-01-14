classdef SMLMModelFitter<interfaces.SEEvaluationProcessor
    % This is a plugin in development. Public has no access to the 
    % run_ functions called in this plugin. For internal user, you need
    % NPC3D to run this plugin.
    methods
        function obj=SMLMModelFitter(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        
        function out=run(obj,p)
            out = runSMLMModelFitter(obj,p);
        end
        function pard=guidef(obj)
            %init
            if ~isdeployed
                addpath(['..' filesep 'ries-private' filesep 'SMLMModelFitter'], ['..' filesep 'ries-private' filesep 'SMLMModelFitter' filesep 'external'])
            end
            pard=SMLMModelFitter_guidef(obj);
        end
    end
end