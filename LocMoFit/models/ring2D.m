classdef ring2D<geometricModel
    % :class:`ring2D` is a 2D model that describes a ring geometry.
    %
    % Geometric parameters:
    %   * `radius`: (nm) the ring radius.
    %
    % Relavent biological structure:
    %   * Top-view projections of the nuclear pore complex.
    %
    % Preview:
	% 	.. note::
	% 		It will be available soon.
	%
	% ..
    %   .. image:: ./images/models/ring2D.PNG
    %       :width: 400
    %   Scale bar: 50 nm.
    methods
        function obj = ring2D(varargin)
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
            obj.dimension = 2;
            obj.listed = true;
            
        end
        
        function [model, p]= reference(obj, par, dx)
        % For details, see :meth:`reference`.
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
        model.channel = ones(size(model.x));
        model.n = ones(size(model.x));

        p = [];
        end
        function derivedPars = getDerivedPars(obj, pars)
            derivedPars = [];
        end
    end
end
