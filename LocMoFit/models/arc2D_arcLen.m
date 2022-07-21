classdef arc2D_arcLen<geometricModel
    % Last update:
    %   21.07.2022
    methods
        function obj = arc2D_arcLen(varargin)
            obj@geometricModel(varargin{:});
            % Define parameters that can be altered during fitting here:
            obj.name = {'arcLength', 'theta'}; % parameter names
            obj.fix = [0 0] ;                                                       % fix to a constant or not
            obj.value = [10 30];                                                    % initial guess
            obj.lb = [-inf -inf];                                                   % relative lower bound
            obj.ub = [inf inf];                                                     % relative upper bound
            obj.min = [5 5];                                                        % absolute lower bound
            obj.max = [30 180];                                                     % absolute upper bound
                       
            % Define other properties here:
            obj.modelType = 'discretized';
            obj.modelTypeOption = {'discretized','continuous'};
            obj.dimension = 2;
            obj.listed = true;
        end
        
        function [model, p]= reference(obj, par, dx)
        % This funtion samples coordinates of the model as reference.
        %
        % Usage:
        %   [model, p]= reference(obj, par, dx)
        %
        % Args:
        %   * 'obj': an object of subclass of :class:`geometricModel`.
        %   * 'par': a structure object. Its fieldnames should be the names
        %   of parameters, and their correspoinding content should be the
        %   parameter values.
        %   * 'dx': sampling rate.
        %
        % Returns:
        %   * 'model': a structure object. Its fieldnames should be x, y,
        %   z, and n, indicating the xyz position amplitude n of the
        %   sampled model points.
        %   * 'p': additional information of the model.

        arcLen = par.arcLength;
        twiceTheta = deg2rad(2*par.theta);
        r = arcLen./twiceTheta;
        % 
        if isempty(obj.ParentObject)
            locsPrecFactor = 1;
        else
            locsPrecFactor = min(obj.ParentObject.locsPrecFactor,5);
        end
%         arcLen = 2*pi*r*theta/360;
        nSample = max(round(arcLen/(dx*locsPrecFactor)),1);
        
        halfRange = twiceTheta./2;
        centroidx = r.*sin(halfRange)/halfRange;
        ang_samplePoints = linspace(-halfRange, halfRange, nSample)';
        model.x = r.*cos(ang_samplePoints)-centroidx;
        model.y = r.*sin(ang_samplePoints);
        [model.x, model.y] = rotcoord2(model.x, model.y, pi/2);
        model.n = ones(size(model.x));

        p.cornerRange = [];
        end
        function derivedPars = getDerivedPars(obj, pars)
            derivedPars = [];
        end
    end
end
