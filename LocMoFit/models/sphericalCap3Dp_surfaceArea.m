classdef sphericalCap3Dp_surfaceArea<parametricModel
    % :class:`sphericalCap3Dp_surfaceArea` describes a geometry of spherical cap in 3D. It describe the same geometry with the same parameterization as :class:`sphericalCap3D_surfaceArea<models.sphericalCap3D_surfaceArea>` but in a parametric form.
	%
	% Geometric parameters:
    %   * `surfaceArea`: (10\ :sup:`4` nm\ :sup:`2`) the surface area of the spherical cap.
    %   * `closeAngle`: (°) the angle from the pole to the edge of the cap.
	%
    % Relavent biological structure:
    %   * mammalian endocytic coat
	%
	% See also:
    %   :class:`sphericalCap3D_surfaceArea<models.sphericalCap3D_surfaceArea>`
	%
    % Preview:
	% 	See :class:`sphericalCap3D_surfaceArea<models.sphericalCap3D_surfaceArea>`
	
    methods
        function obj = sphericalCap3Dp_surfaceArea(varargin)
            obj@parametricModel(varargin{:});
            % Define the default argument values here in the constructor.
            obj.name = {'surfaceArea', 'closeAngle'};
            obj.fix = [0 0] ;
            obj.value = [5 0];
            obj.lb = [-inf 0];
            obj.ub = [inf 90];
            obj.min = [-50 0];
            obj.max = [50 180];
            
            % Define other properties here:
            obj.modelType = 'discretized';
            obj.modelTypeOption = {'discretized', 'continuous'};
            obj.dimension = 3;
            obj.listed = true;
            
            % Specific to a parametric model
        end

        function [model, p]= definedModel(obj, u, v, par, dx)
            [U,V] = meshgrid(u,v);
			%% Get parameters
            U = U(:); V = V(:);
            if isempty(obj.ParentObject)||isempty(obj.ParentObject.locsPrecFactor)
                    locsPrecFactor = 1;
                else
                    locsPrecFactor = obj.ParentObject.locsPrecFactor;
            end

            if par.r ~= inf
                allC = sqrt(par.r.^2-U.^2);
                x = allC.*cos(V);
                y = allC.*sin(V);
                z = U;
                p.zOffset = par.r*0.5*(sin(par.closeAngle+pi/2)+1);
            else
                x = U.*cos(V);
                y = U.*sin(V);
                z = zeros(size(U));
                p.zOffset = 0;
            end
            AllPoints = [x y z];
            %             AllPoints = unique(AllPoints, 'rows');
            [~,d] = rangesearch(AllPoints,AllPoints,locsPrecFactor*1);
            density = cellfun(@(x)length(x), d);
            n = 1./density;

            %% Get the final model
            % Scale the unit sphere to radius
            model.x = AllPoints(:,1);
            model.y = AllPoints(:,2);
            model.z = AllPoints(:,3);      % upside-down the coordinates

            % move the center from the sphere origin to the centroid of mass
            %             p.zOffset = median(z);
            
            model.z = model.z - p.zOffset;
            model.n = n;

            % validate the distance
            %             meanD = validation_density(model.x, model.y, model.z);
        end
        
        function [u,v] = getParVector(obj,par,dx)
            if isempty(obj.ParentObject)||isempty(obj.ParentObject.locsPrecFactor)            %             rimPos = par.r.*sin(-par.closeAngle);
                locsPrecFactor = 1;
            else
                locsPrecFactor = obj.ParentObject.locsPrecFactor;
            end
            minD = dx*locsPrecFactor;
            if par.r == inf
                uStep = ceil(par.rc_max./minD);
                u = linspace(minD, par.rc_max, uStep-1);
                maxC = 2.*pi.*par.rc_max;
            else
                rim = abs(par.r.*sin(par.closeAngle));
                arcLen = pi*par.r*(par.closeAngle/pi); % 2pi*R(ratio) devided by 2

                uStep = abs(ceil(arcLen./minD));
                % in polar coordinates
                dTheta = linspace(pi/2, pi/2-par.closeAngle, uStep);
                dTheta(1) = [];
                u = par.r.*sin(dTheta);

                if par.closeAngle>=pi/2
                    maxC = abs(2*pi*par.r);
                else
                    maxC = 2*pi*rim;
                end
            end
            vStep = ceil(maxC./minD);
            v = deg2rad(linspace(0, 360, vStep));

        end
        
        
        function par = convertPar(obj, par)
            par_.surfaceArea = par.surfaceArea*1e+4;
            % Get the close angle from the input
            par_.closeAngle = par.closeAngle;
            lNeg = par_.closeAngle<0;

            par_.closeAngle = deg2rad(abs(par_.closeAngle));
            
            par_.r = sqrt(abs(par_.surfaceArea)/(2*pi*(1-cos(par_.closeAngle))));         
            par_.r(lNeg) = -par_.r(lNeg);
            
            if abs(par_.r) == inf
                par_.rc_max = sqrt(par_.surfaceArea/pi);
            end
            par = par_;
        end
        
        function derivedPars = getDerivedPars(obj, pars)
            % For details, see :meth:`getDerivedPars`.
			
            derivedPars.realSurfaceArea = pars.surfaceArea.*1e+4;
            derivedPars.curvature = 1./sqrt(derivedPars.realSurfaceArea./(2.*pi.*(1-cos(deg2rad(pars.closeAngle)))));
            if pars.closeAngle < 0
                derivedPars.curvature = -derivedPars.curvature;
                derivedPars.realCloseAngle = abs(pars.closeAngle);
            else
                derivedPars.realCloseAngle = pars.closeAngle;
            end
            derivedPars.radius = 1/derivedPars.curvature;
            
            if derivedPars.realCloseAngle <= 90
                derivedPars.projectionArea = pi*(derivedPars.radius*sin(deg2rad(derivedPars.realCloseAngle)))^2;
            else
                derivedPars.projectionArea = pi*derivedPars.radius^2;
            end
            derivedPars.coverageFraction = (1-cos(deg2rad(derivedPars.realCloseAngle)))./2;
            
            % surface area (one side free)
            signRadius = sign(derivedPars.radius);
            radiusOneSideFree = derivedPars.radius-signRadius*pars.variation;
            derivedPars.areaOneSideFree = radiusOneSideFree.^2*(2.*pi.*(1-cos(deg2rad(pars.closeAngle))));
            
            %% Final version used in publication
            derivedPars.openRadius = abs(derivedPars.radius*sin(deg2rad(180-derivedPars.realCloseAngle)));
            if ~isempty(obj.ParentObject)&&~isempty(obj.ParentObject.ParentObject)
                locMoFitter = obj.ParentObject.ParentObject;
                modID = obj.ParentObject.ID;
                derivedPars.basePos = derivedPars.radius*0.5*(sin(deg2rad(derivedPars.realCloseAngle)+pi/2)+1);
            end
        end
    end
end