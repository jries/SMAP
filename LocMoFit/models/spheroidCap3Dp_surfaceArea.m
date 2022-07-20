classdef spheroidCap3Dp_surfaceArea<parametricModel
    % Describing endocytic coat proteins as molecules covering a part of
    % sphere with an angle indicating the closed part.
    %
    % Last update:
    %   18.07.2021
    methods
        function obj = spheroidCap3Dp_surfaceArea(varargin)
            obj@parametricModel(varargin{:});
            % Define the default argument values here in the constructor.
            obj.name = {'surfaceArea', 'closeAngle', 'flattening'};
            obj.fix = [0 0 0] ;
            obj.value = [5 0 0];
            obj.lb = [-inf 0 -inf];
            obj.ub = [inf 90 inf];
            obj.min = [0 -90 0];
            obj.max = [50 90 0.7];
            
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
            if isempty(obj.ParentObject)||isempty(obj.ParentObject.locsPrecFactor)
                    locsPrecFactor = 1;
                else
                    locsPrecFactor = obj.ParentObject.locsPrecFactor;
            end

            if par.r ~= inf
                allA = sqrt(par.A.^2-U.^2);
%                 allA;
                allC = allA.*(1-par.flattening);
                x = allA.*cos(V);
                y = allC.*sin(V);
                z = U;
                p.zOffset = par.r*0.5*(sin(par.closeAngle+pi/2)+1);
            else
                %!!! yet to be modified
                avgR = U;
                c = 2.*avgR*(1-par.flattening)./(2-par.flattening);
                a = 2.*avgR-c;
                x = a.*cos(V);
                y = c.*sin(V);
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
            %             rimPos = par.r.*sin(-par.closeAngle);
            if isempty(obj.ParentObject)||isempty(obj.ParentObject.locsPrecFactor)
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
%                 C = par.C;
                A = par.A;

                rim = abs(A.*sin(par.closeAngle));
                arcLen = pi*A*(par.closeAngle/pi); % 2pi*R(ratio) devided by 2

                uStep = abs(ceil(arcLen./minD));
                % in polar coordinates
                dTheta = linspace(pi/2, pi/2-par.closeAngle, uStep);
                dTheta(1) = [];
                u = A.*sin(dTheta);

                if par.closeAngle>=pi/2
                    maxC = abs(2*pi*A);
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
            par_.flattening = par.flattening;
            lNeg = par_.closeAngle<0;

            par_.closeAngle = deg2rad(abs(par_.closeAngle));
            
            par_.r = sqrt(abs(par_.surfaceArea)/(2*pi*(1-cos(par_.closeAngle))));         
            par_.r(lNeg) = -par_.r(lNeg);

            par_.C = 2.*par_.r*(1-par.flattening)./(2-par.flattening);
            par_.A = 2.*par_.r-par_.C;
            
            if abs(par_.r) == inf
                par_.rc_max = sqrt(par_.surfaceArea/pi);
                par_.Cc_max = 2.*par_.rc_max*(1-par.flattening)./(2-par.flattening);
                par_.Ac_max = 2.*par_.r-par_.C;
            end
            par = par_;
        end
        
        function derivedPars = getDerivedPars(obj, pars)
            % Get extra parameters derived from the fit parameters.
            %
            % Exports a empty variable when no derived parameters.
            % Last update:
            %   08.11.2021
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