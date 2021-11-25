classdef CME3DSphereCoverageCur_discrete<geometricModel
    % Describing endocytic coat proteins as molecules covering a part of
    % sphere with an angle indicating the closed part.
    methods
        function obj = CME3DSphereCoverageCur_discrete(varargin)
            obj@geometricModel(varargin{:});
            % Define the default argument values here in the constructor.
            obj.name = {'curvature', 'closeAngle'};
            obj.fix = [0 0] ;
            obj.value = [0 0];
            obj.lb = [-inf 0];
            obj.ub = [inf 90];
            obj.min = [5e-1 -90];
            obj.max = [30 90];
            obj.modelType = 'discrete';
            obj.dimension = 3;
        end

        function [model, p]= reference(obj, par, dx)
            %% Get parameters
            % Get radius based on the input 1000X curvature
            cur = par.curvature;
            r = (cur./1000).^-1; % the 1000 here is a factor to adjust the range. Originally it was too small.
            
            % Get the close angle from the input
            closeAngle = deg2rad(par.closeAngle);
            transition_cur = 0.005;
            if abs(r) > transition_cur^-1
                r_area = transition_cur^-1;
                % area_cap = 2*pi*r*h; h = area_cap/2*pi*r;
                % r*sin(theta)=r-h; h = -r(sin(theta)-1);
                % area_cap = 2*pi*r*(-r(sin(theta)-1)) = -2*pi*r^2*(sin(theta)-1)
                area = 2*pi*r_area^2*(1-cos(closeAngle+pi/2));
                closeAngle = acos(-area/(2*pi*r^2)+1)-pi/2;
            end
            if r == inf
                rc_max = sqrt(area/pi);
            end
            %% Evenly distribute points using a Fibonacci sphere
            % the sampling number should be propotional to surface area
            
            samplingFactor = round(15/dx)*obj.ParentObject.locsPrecFactor/10;      % change this if you need more points
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
            
            
            if r < inf
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
            p.closeAngle_real = rad2deg(closeAngle);
            model.z = model.z - p.zOffset;
            model.n = ones(size(model.x));
            
            meanD = validation_density(model.x, model.y, model.z);
        end
    end
end
function meanD = validation_density(x,y,z)
    qPoints = [x y z];
    [~,d] = knnsearch(qPoints, qPoints, 'k', 2);
    meanD = mean(d(:,2));
end