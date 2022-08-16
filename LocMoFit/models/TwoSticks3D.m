classdef TwoSticks3D<geometricModel
    methods
        function obj = TwoSticks3D(varargin)
            obj@geometricModel(varargin{:});
            % Define parameters that can be altered during fitting here:
            obj.name = {'Tlength', 'distance','distortion'}; % parameter names
            obj.fix = [0 0 1] ;                                                   % fix to a constant or not
            obj.value = [80 110 35];                                                % initial guess
            obj.lb = [-10 -10 0];                                                 % relative lower bound
            obj.ub = [10 10 60];                                                   % relative upper bound
            obj.min = [5 5 0];                                                    % absolute lower bound
            obj.max = [70 50 50];                                                  % absolute upper bound
            
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
        
        %circumference = 2*pi*par.distance;
        
        if isempty(obj.ParentObject.locsPrecFactor)
            locsPrecFactor = 1;
        else
            locsPrecFactor = obj.ParentObject.locsPrecFactor;
        end
        minDist = locsPrecFactor*dx;      % change this if you need more points

        % corner:
        theta = [0 pi];
        %theta(end) = [];
        
        % assign the distance for each copy
        rho = repelem(par.distance', size(theta,1));
        
        % convert from polar coordinates
        [oneRingX, oneRingY]= pol2cart(theta(:), rho(:));
        
        %pointPerRing = length(oneRingX);
        allRingZ = linspace(-par.Tlength/2, par.Tlength/2, par.Tlength/minDist);
        allRingX = repmat([oneRingX; oneRingX; oneRingX], [length(allRingZ) 1]);
        allRingY = repmat([oneRingY-par.distortion;oneRingY;oneRingY+par.distortion], [length(allRingZ) 1]);
        
        allRingZ = repelem(allRingZ', length([oneRingX;oneRingX;oneRingX]));

        model.x = allRingZ;
        model.y = allRingX;
        model.z = allRingY;

        model.channel = ones(size(model.x));
        model.n = ones(size(model.x));
        p = [];
        end
        function derivedPars = getDerivedPars(obj, pars)
            derivedPars = [];
        end
    end
end
