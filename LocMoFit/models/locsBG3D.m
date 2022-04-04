classdef locsBG3D<geometricModel
    
    % log
    %   - 201229: change the sign of the ring twist
    methods
        function obj = locsBG3D(varargin)
            obj@geometricModel(varargin{:});
            % Define parameters that can be altered during fitting here:
            obj.name = {'density', 'depth'};                               % parameter names
            obj.fix = [0 1];                                               % fix to a constant or not
            obj.value = [20 200];                                          % initial guess
            obj.lb = [-inf -inf];                                          % relative lower bound
            obj.ub = [inf inf];                                            % relative upper bound
            obj.min = [0 100];                                             % absolute lower bound
            obj.max = [inf 1000];                                          % absolute upper bound
                       
            % Define other properties here:
            obj.modelType = 'discretized';
            obj.modelTypeOption = {'discretized','continuous'};
            obj.dimension = 3;
        end
        
        function [model, p]= reference(obj, par, dx)
            
        % Sample coordinates of the model as reference.
        % --- Syntax ---
        % [model, p]= reference(obj, par, dx)
        % --- Arguments ---
        % -- Input --
        % obj:
        % par: a structure object. Its fieldnames should be the names of
        % parameters, and their correspoinding content should be the
        % parameter values.
        % dx: sampling rate.
        % -- Output --
        % model: a structure object. Its fieldnames should be x, y, z, and
        % n, indicating the xyz position amplitude n of the sampled model
        % points.
        % p: additional information of the model.
        
        roiSize = obj.ParentObject.ParentObject.roiSize;
        copyNumber = par.density.*roiSize./1000;
        points = (rand([copyNumber, 3])-0.5) .* [roiSize roiSize par.depth];
        
        model.x = points(:,1);
        model.y = points(:,2);
        model.z = points(:,3);
        model.n = ones(size(model.x));
        %     likelihood = gaussDist(model,double(locs.ynmrot), double(locs.xnmrot), double(locs.znm), 8,8,10);
        %     likelihoodBoun = gaussDist(model,0,0,0, 16,16,20);

        p = [];
        end
        function derivedPars = getDerivedPars(obj, pars)
            derivedPars = [];
        end
    end
end
