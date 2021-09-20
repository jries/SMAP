classdef cylinder3D<geometricModel
    methods
        function obj = cylinder3D(varargin)
            obj@geometricModel(varargin{:});
            % Define parameters that can be altered during fitting here:
            obj.name = {'height', 'radius'}; % parameter names
            obj.fix = [0 0] ;                                                   % fix to a constant or not
            obj.value = [30 15];                                                % initial guess
            obj.lb = [-10 -10];                                                 % relative lower bound
            obj.ub = [10 10];                                                   % relative upper bound
            obj.min = [5 5];                                                    % absolute lower bound
            obj.max = [70 50];                                                  % absolute upper bound
            
            % Define discrite parameters here:
            
            % Define other properties here:
            obj.modelType = 'discretized';
            obj.modelTypeOption = {'discretized', 'continous'};
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
        
        circumference = 2*pi*par.radius;
        
        if isempty(obj.ParentObject.locsPrecFactor)
            locsPrecFactor = 1;
        else
            locsPrecFactor = obj.ParentObject.locsPrecFactor;
        end
        minDist = locsPrecFactor*dx;      % change this if you need more points

        % corner:
        theta = linspace(0,2*pi,ceil(circumference/minDist));
        theta(end) = [];
        
        % assign the radius for each copy
        rho = repelem(par.radius', size(theta,1));
        
        % convert from polar coordinates
        [oneRingX, oneRingY]= pol2cart(theta(:), rho(:));
        pointPerRing = length(oneRingX);
        allRingZ = linspace(-par.height/2, par.height/2, par.height/minDist);
        
        allRingX = repmat(oneRingX, [length(allRingZ) 1]);
        allRingY = repmat(oneRingY, [length(allRingZ) 1]);
        
        allRingZ = repelem(allRingZ', pointPerRing);

        model.x = allRingX;
        model.y = allRingY;
        model.z = allRingZ;

        model.channel = ones(size(model.x));
        model.n = ones(size(model.x));
        p = [];
        end
        function derivedPars = getDerivedPars(obj, pars)
            derivedPars = [];
        end
    end
end
