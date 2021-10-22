classdef parametricModel<geometricModel
    % :class:`parametricModel` is a subclass of `geometricModel`. This class
    % allows the user to define their model in a parametric manner.
    %
    methods (Sealed = true)
        function obj = parametricModel(varargin)
            obj@geometricModel(varargin{:})
        end
        function [model,p] = reference(obj, par, dx)
            % Please do not modify this part. If you have any requst
            % regarding this part, please contact us.
            
            %% The internal conversion of parameters
            par = obj.convertPar(par);
            
            %% Get parametric vectors
            
            % Depeding on the dimensionality of the model
            switch obj.dimension
                case 2
                    u = obj.getParVector(par, dx);
                case 3
                    [u,v] = obj.getParVector(par, dx);
            end
            
            %% Sample the model according to the parametric vectors
            [model, p] = obj.definedModel(u,v,par,dx);
        end
    end
    methods (Abstract = true)
        % Define your geometric model here.
        [model, p] = definedModel(obj,varargin)
        
        % Define how to create the parametric vectors u and v here.
        [u,v] = getParVector(obj,varargin)
    end
    
    methods
        function par = convertPar(obj, par)
            % Define the internal conversion of parameters. Skip this part
            % if no conversion is necessary.
            par = par;
        end
    end
end
