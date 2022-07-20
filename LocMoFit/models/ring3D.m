classdef ring3D<geometricModel
    % :class:`dualRingModel` is the dual ring model used in the LocMoFit
    % manuscript for describing Nup96-labeled NPCs.
    %
    % Log:
    %   201229: change the sign of the ring twist
    methods
        function obj = ring3D(varargin)
            obj@geometricModel(varargin{:});
            % Define parameters that can be altered during fitting here:
            obj.name = {'radius'}; % parameter names
            obj.fix = 1;           % fix to a constant or not
            obj.value = 53.7;      % initial guess
            obj.lb = 0;            % relative lower bound
            obj.ub = 0;            % relative upper bound
            obj.min = 0;           % absolute lower bound
            obj.max = 100;         % absolute upper bound
            
            % Define other properties here:
            obj.modelType = 'discretized';
            obj.modelTypeOption = {'discretized','continuous'};
            obj.dimension = 3;
            obj.listed = true;
            
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
        
                
        if isempty(obj.ParentObject.locsPrecFactor)
                locsPrecFactor = 1;
        else
            locsPrecFactor = obj.ParentObject.locsPrecFactor;
        end
        
        minD = locsPrecFactor*dx;
        
        % corner:
        circumference = 2.*pi.*par.radius;
        cornerNum = circumference./minD;
        cornerPos = linspace(0,2*pi,cornerNum+1);
        cornerPos = cornerPos(1:end-1);

        [allCopiesX,allCopiesY] = pol2cart(cornerPos(:), par.radius(:));
        
        model.x = allCopiesX;
        model.y = allCopiesY;
        model.z = zeros(size(model.x));
        model.channel = ones(size(model.x));
        model.n = ones(size(model.x));

        p = [];
        end
        function derivedPars = getDerivedPars(obj, pars)
            derivedPars = [];
        end
    end
end
