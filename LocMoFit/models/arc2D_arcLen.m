classdef arc2D_arcLen<geometricModel
    % :class:`arc2D_arcLen` is a 2D model that describes an arc geometry.
    % It describes the same geometry as by :class:`arc2D<models.arc2D>` but
    % with a different parameterization.
    %
    % Geometric parameters:
    %   * `arcLength`: (nm) the length of the arc.
    %   * `theta`: (degree) the closing anlge of the arc.
    %
    % Relavent biological structure:
    %   * Cross-section of a clathrin coat
    %
    % See also:
    %   :class:`arc2D<models.arc2D>`
    
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
        % For details, see :meth:`reference`.

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
