classdef bucket2D<geometricModel
    % :class:`bucket2D` is a 2D model that describes a bucket geometry. It
    % creates a bucket based on an arc. The bucket is created to contain
    % the arc.
    %
    % Geometric parameters:
    %   * `radius`: (nm) the radius of the ring where the arc is derived.
    %   * `theta`: (°) the closing anlge of the arc.
    %
    % See also:
    %   :class:`arc2D<models.arc2D>`
    methods
        function obj = bucket2D(varargin)
            obj@geometricModel(varargin{:});
            % Define parameters that can be altered during fitting here:
            obj.name = {'radius', 'theta'}; % parameter names
            obj.fix = [0 0] ;                                                       % fix to a constant or not
            obj.value = [10 60];                                                    % initial guess
            obj.lb = [-inf -inf];                                                   % relative lower bound
            obj.ub = [inf inf];                                                     % relative upper bound
            obj.min = [5 5];                                                        % absolute lower bound
            obj.max = [30 360];                                                     % absolute upper bound
                       
            % Define other properties here:
            obj.modelType = 'discretized';
            obj.modelTypeOption = {'discretized','continuous'};
            obj.dimension = 2;
            obj.listed = true;
        end
        
        function [model, p]= reference(obj, par, dx)
        % For details, see :meth:`reference`.
        
        % parameters
        r = par.radius;
        theta = par.theta;

        % 
        if isempty(obj.ParentObject)
            locsPrecFactor = 1;
        else
            locsPrecFactor = min(obj.ParentObject.locsPrecFactor,5);
        end
        
        halfTheta = deg2rad(theta/2);
        xp = r*cos(halfTheta);
        yp = r*sin(halfTheta);
        if theta<=180
            xAnchor = [xp, r, r, xp];
            yAnchor = [yp, yp, -yp, -yp];
        else
            xAnchor = [xp, xp, r, r, xp, xp];
            yAnchor = [yp, r, r, -r, -r, -yp];
        end
        arcLen = arclength(xAnchor,yAnchor);
        nSample = max(round(arcLen/(dx*locsPrecFactor)),1);
        
        pt = interparc(nSample,xAnchor,yAnchor,'linear');
        
        centroidx = mean(pt(:,1));
        model.x = pt(:,1)-centroidx;
        model.y = pt(:,2);
        model.n = ones(size(model.x));
        
        p = [];
        end
        function derivedPars = getDerivedPars(obj, pars)
            % For details, see :meth:`getDerivedPars`.
            derivedPars = [];
        end
    end
end
