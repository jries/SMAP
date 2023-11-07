classdef ellipse2D<geometricModel
    % :class:`ellipse2D` describes a ellipse geometry in 2D.
    %
    % Geometric parameters:
    %   * `a`: (nm) the axis along the x-axis.
	%   * `b`: (nm) the axis along the y-axis.
    %
    % Preview:
	% 	.. note::
	% 		It will be available soon.
	%
	% ..
    %   .. image:: ./images/models/ellipse2D.PNG
    %       :width: 400
    %   Scale bar: 50 nm.
    methods
        function obj = ellipse2D(varargin)
            obj@geometricModel(varargin{:});
            % Define parameters that can be altered during fitting here:
            obj.name = {'a','b'}; % parameter names
            obj.fix = [1 1];           % fix to a constant or not
            obj.value = [40 50];      % initial guess
            obj.lb = [-inf -inf];            % relative lower bound
            obj.ub = [inf inf];            % relative upper bound
            obj.min = [0 0];           % absolute lower bound
            obj.max = [100 100];         % absolute upper bound
            
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
        h = ((par.a-par.b)./(par.a+par.b)).^2;
        circumference = pi.*(par.a+par.b).*(1+3.*h./(10+sqrt(4-3*h)));
        cornerNum = circumference./minD;
		
		theta = linspace(0,2*pi,cornerNum+1);
        theta = theta(1:end-1);
        
        oneRingX = par.a*cos(theta);
        oneRingY= par.b*sin(theta);
        
        model.x = oneRingX.';
        model.y = oneRingY.';

        model.channel = ones(size(model.x));
        model.n = ones(size(model.x));

        p = [];
        end
        function derivedPars = getDerivedPars(obj, pars)
            derivedPars = [];
        end
    end
end
