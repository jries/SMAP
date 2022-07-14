classdef sphericalCap3D_surfaceArea_par<parametricModel
    % Describing endocytic coat proteins as molecules covering a part of
    % sphere with an angle indicating the closed part.
    %
    % Last update:
    %   14.07.2021
    methods
        function obj = sphericalCap3D_surfaceArea_par(varargin)
            obj@parametricModel(varargin{:});
            % Define the default argument values here in the constructor.
            obj.name = {'surfaceArea', 'closeAngle'};
            obj.fix = [0 0] ;
            obj.value = [0.5e+5 0];
            obj.lb = [-inf 0];
            obj.ub = [inf 90];
            obj.min = [-5e+5 -90];
            obj.max = [5e+5 90];
            
            % Define other properties here:
            obj.modelType = 'discretized';
            obj.modelTypeOption = {'discretized', 'continuous'};
            obj.dimension = 3;
            obj.listed = true;
            
            % Specific to a parametric model
        end

        function [model, p]= definedModel(obj, u, v, par, dx)
            %% Get parameters
            [U,V] = meshgrid(u,v);
            U = U(:); V = V(:);
            allC = sqrt(par.r.^2-U.^2);
            x = allC.*cos(V);
            y = allC.*sin(V);
            z = U;
            
            AllPoints = [x y z];
%             AllPoints = unique(AllPoints, 'rows');
            
            if isempty(obj.ParentObject)||isempty(obj.ParentObject.locsPrecFactor)
                locsPrecFactor = 1;
            else
                locsPrecFactor = obj.ParentObject.locsPrecFactor;
            end
            [~,d] = rangesearch(AllPoints,AllPoints,locsPrecFactor*1.5);
            density = cellfun(@(x)length(x), d);
            n = 1./density;
            %% Get the final model
            % Scale the unit sphere to radius
            model.x = AllPoints(:,1);
            model.y = AllPoints(:,2);
            model.z = AllPoints(:,3);      % upside-down the coordinates
            
            % move the center from the sphere origin to the centroid of mass 
%             p.zOffset = median(z);
             p.zOffset = par.r*0.5*(sin(-par.closeAngle)+1);
            model.z = model.z - p.zOffset;
            model.n = n;
            
            % validate the distance
%             meanD = validation_density(model.x, model.y, model.z);
        end
        
        function [u,v] = getParVector(obj,par,dx)
%             rimPos = par.r.*sin(-par.closeAngle);
            rim = par.r.*cos(-par.closeAngle);
            arcLen = 2*pi*par.r*((0.5*pi-par.closeAngle)./pi);
            
            
            if isempty(obj.ParentObject)||isempty(obj.ParentObject.locsPrecFactor)
                locsPrecFactor = 1;
            else
                locsPrecFactor = obj.ParentObject.locsPrecFactor;
            end
            
            uStep = ceil(arcLen./(dx*locsPrecFactor));
            dTheta = linspace(pi/2-0.05, -par.closeAngle, uStep);
            u = par.r.*sin(dTheta);
            
            if par.closeAngle>=0
                maxC = 2*pi*par.r;
            else
                maxC = 2*pi*rim;
            end
            
            vStep = ceil(maxC./(dx*locsPrecFactor));
            v = deg2rad(linspace(0, 360, vStep));
        end
        
        
        function par = convertPar(obj, par)
            par_.surfaceArea = par.surfaceArea*1e+4;
            % Get the close angle from the input
            par_.closeAngle = par.closeAngle;
            lNeg = par_.closeAngle<-90;
            if lNeg
                par_.closeAngle = abs(par_.closeAngle)-180;
            end
            par_.closeAngle = deg2rad(par_.closeAngle);
            
            par_.r = sqrt(abs(par_.surfaceArea)/(2*pi*(1-cos(par_.closeAngle+pi/2))));         
            par_.r(lNeg) = -par_.r(lNeg);
            
            if abs(par_.r) == inf
                par_.rc_max = sqrt(par_.surfaceArea/pi);
            end
            par = par_;
        end
        
        function derivedPars = getDerivedPars(obj, pars)
            % Exports a empty variable when no derived parameters.
            derivedPars.realSurfaceArea = pars.surfaceArea.*1e+4;
            derivedPars.curvature = 1./sqrt(derivedPars.realSurfaceArea./(2.*pi.*(1-cos(deg2rad(90+pars.closeAngle)))));
            if pars.closeAngle < -90
                derivedPars.curvature = -derivedPars.curvature;
                derivedPars.realCloseAngle = abs(pars.closeAngle)-180;
            else
                derivedPars.realCloseAngle = pars.closeAngle;
            end
            derivedPars.radius = 1/derivedPars.curvature;
            if derivedPars.realCloseAngle <= 0
                derivedPars.projectionArea = pi*(derivedPars.radius*cos(deg2rad(-derivedPars.realCloseAngle)))^2;
            else
                derivedPars.projectionArea = pi*derivedPars.radius^2;
            end
            derivedPars.coverageFraction = (1-cos(deg2rad(90+derivedPars.realCloseAngle)))./2;
            
            % surface area (one side free)
            signRadius = sign(derivedPars.radius);
            radiusOneSideFree = derivedPars.radius-signRadius*pars.variation;
            derivedPars.areaOneSideFree = radiusOneSideFree.^2*(2.*pi.*(1-cos(deg2rad(90+pars.closeAngle))));
        end
    end
end