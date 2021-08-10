classdef CME3DSphereCoverageArea_discrete<geometricModel
    % Describing endocytic coat proteins as molecules covering a part of
    % sphere with an angle indicating the closed part.
    methods
        function obj = CME3DSphereCoverageArea_discrete(varargin)
            obj@geometricModel(varargin{:});
            % Define the default argument values here in the constructor.
            obj.name = {'surfaceArea', 'closeAngle'};
            obj.fix = [0 0] ;
            obj.value = [0.5e+5 0];
            obj.lb = [-inf 0];
            obj.ub = [inf 90];
            obj.min = [-5e+5 -90];
            obj.max = [5e+5 90];
            
            % Define other properties here:
            obj.modelType = 'discrete';
            obj.modelTypeOption = {'discretized', 'continuous'};
            obj.dimension = 3;
        end

        function [model, p]= reference(obj, par, dx)
            %% Get parameters
            surfaceArea = par.surfaceArea*1e+4;
            % Get the close angle from the input
            closeAngle = par.closeAngle;
            lNeg = closeAngle<-90;
            if lNeg
                closeAngle = abs(closeAngle)-180;
            end
            closeAngle = deg2rad(closeAngle);
            
            r = sqrt(abs(surfaceArea)/(2*pi*(1-cos(closeAngle+pi/2))));         
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
            samplesUb = round((sin(closeAngle)+1)./offset-0.5);
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

function meanD = validation_density(x,y,z)
    qPoints = [x y z];
    [~,d] = knnsearch(qPoints, qPoints, 'k', 2);
    meanD = mean(d(:,2));
end