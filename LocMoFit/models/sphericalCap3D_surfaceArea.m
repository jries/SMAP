classdef sphericalCap3D_surfaceArea<geometricModel
    % :class:`sphericalCap3D_surfaceArea` describes the geometry of a spherical cap in 3D.
	%
	% Geometric parameters:
    %   * `surfaceArea`: (10\ :sup:`4` nm\ :sup:`2`) the surface area of the spherical cap.
    %   * `closeAngle`: (°) the angle from the pole to the edge of the cap.
	%
    % Relavent biological structure:
    %   * mammalian endocytic coat
	%
    % Preview:
    %   .. image:: ./images/models/sphericalCap3D_surfaceArea.png
    %       :width: 400
    %   Scale bar: 50 nm.

    methods
        function obj = sphericalCap3D_surfaceArea(varargin)
            obj@geometricModel(varargin{:});
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
        end

        function [model, p]= reference(obj, par, dx)
		% For details, see :meth:`reference`.
		
            %% Get parameters
            surfaceArea = par.surfaceArea*1e+4;
            % Get the close angle from the input
            closeAngle = par.closeAngle;
            lNeg = closeAngle<0;
            if lNeg
                closeAngle = abs(closeAngle);
            end
            closeAngle = deg2rad(closeAngle);
            
            r = sqrt(abs(surfaceArea)/(2*pi*(1-cos(closeAngle))));         
            r(lNeg) = -r(lNeg);
            
            if abs(r) == inf
                rc_max = sqrt(surfaceArea/pi);
            end
            %% Evenly distribute points using a Fibonacci sphere
            % the sampling number should be propotional to surface area
%            aaa(end+1) = ;
            if isempty(obj.ParentObject.locsPrecFactor)
                locsPrecFactor = 1;
            else
                locsPrecFactor = obj.ParentObject.locsPrecFactor;
            end
            samplingFactor = (33.9/(locsPrecFactor*dx))^2;      % change this if you need more points
            if r ~= inf
                samplingScale = round((4*pi*r^2)/(4*pi)*0.01);
            else
                samplingScale = round(pi*rc_max^2/pi*0.01);
            end
            
            if samplingScale == 0
                samplingScale = 1;
            end
            samples = 1 * samplingFactor * samplingScale;       % number of sampled points.
            
            % Unit Fibonacci sphere
            rnd = 1;
            offset = 2./samples;
            gr = (sqrt(5.0) + 1.0) / 2.0;                       % golden ratio = 1.6180339887498948482
            ga = (2.0 - gr) * (2.0*pi);                         % golden angle = 2.39996322972865332
            increment = ga;
            samplesUb = round((sin(closeAngle-pi/2)+1)./offset-0.5);
            k = 1:samplesUb;
            
            
            if abs(r) < inf
                z = (k+0.5).* offset - 1;
                rs = zeros(size(k));
                if samplesUb>samples/2
                    rs(1:(round(samples/2))) = (1 - z(1:(round(samples/2))).^2).^0.5;
                    rs(round(samples/2)+1:samplesUb) = rs(round(samples/2):-1:2*round(samples/2)-samplesUb+1);
                else
                    rs(1:samplesUb) = (1 - z(1:samplesUb).^2).^0.5;          % reduced the computational cost by
                end
                phi = rem((k + rnd), samples) .* increment; %
                
                x = cos(phi) .* rs(k);
                y = sin(phi) .* rs(k);
                
                y = real(y);
            else
                phi = (1:samples).* increment; %
                k = linspace(0,1,samples);
                rc = sqrt(k);
                x = cos(phi) .* rc;
                y = sin(phi) .* rc;
                z = zeros(size(x));
                r = rc_max;
            end
            
            
            %% Get the final model
            % Scale the unit sphere to radius
            model.y = real(r*y)';
            model.x = real(r*x)';
            model.z = -r*z';      % upside-down the coordinates
            
            % move the center from the sphere origin to the centroid of mass 
            p.zOffset = median(model.z);
            model.z = model.z - p.zOffset;
            model.n = ones(size(model.x));
            
            % validate the distance
%             meanD = validation_density(model.x, model.y, model.z);
        end
        
        function derivedPars = getDerivedPars(obj, pars)
            % Get extra parameters derived from the fit parameters.
            %
            % Exports a empty variable when no derived parameters.
            % Last update:
            %   02.07.2023
            derivedPars.realSurfaceArea = pars.surfaceArea.*1e+4;
            derivedPars.curvature = 1./sqrt(derivedPars.realSurfaceArea./(2.*pi.*(1-cos(deg2rad(pars.closeAngle)))));
            if pars.closeAngle < 0
                derivedPars.curvature = -derivedPars.curvature;
                derivedPars.realCloseAngle = abs(pars.closeAngle);
            else
                derivedPars.realCloseAngle = pars.closeAngle;
            end
            derivedPars.radius = 1/derivedPars.curvature;
            
            if derivedPars.realCloseAngle <= 0
                derivedPars.projectionArea = pi*(derivedPars.radius*cos(deg2rad(-derivedPars.realCloseAngle-90)))^2;
            else
                derivedPars.projectionArea = pi*derivedPars.radius^2;
            end
            derivedPars.coverageFraction = (1-cos(deg2rad(derivedPars.realCloseAngle)))./2;
            
            % surface area (one side free)
            signRadius = sign(derivedPars.radius);
            radiusOneSideFree = derivedPars.radius-signRadius*pars.variation;
            derivedPars.areaOneSideFree = radiusOneSideFree.^2*(2.*pi.*(1-cos(deg2rad(pars.closeAngle))));
            
            %% Final version used in publication
            derivedPars.closingAngle_pub = derivedPars.realCloseAngle;
            derivedPars.openRadius = abs(derivedPars.radius*sin(deg2rad(180-derivedPars.closingAngle_pub)));
            if ~isempty(obj.ParentObject)&&~isempty(obj.ParentObject.ParentObject)
                locMoFitter = obj.ParentObject.ParentObject;
                modID = obj.ParentObject.ID;
                derivedPars.basePos = derivedPars.radius*sin(deg2rad(derivedPars.realCloseAngle-90))';
            end
        end
    end
end

function meanD = validation_density(x,y,z)
    qPoints = [x y z];
    [~,d] = knnsearch(qPoints, qPoints, 'k', 2);
    meanD = mean(d(:,2));
end
