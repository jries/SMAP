classdef CME3DSphereCoverageCur<geometricModel
    % Describing endocytic coat proteins as molecules covering a part of
    % sphere with an angle indicating the closed part.
    methods
        function obj = CME3DSphereCoverageCur
            % Define the default argument values here in the constructor.
            obj.name = {'curvature', 'closeAngle'};
            obj.fix = [0 0] ;
            obj.value = [0 0];
            obj.lb = [-inf 0];
            obj.ub = [inf 90];
            obj.min = [5e-1 -90];
            obj.max = [30 90];
            obj.modelType = 'continuous';
            obj.dimension = 3;
        end

        function [model, p]= reference(obj, par, dx)
            
            cur = par.curvature;
            r = (cur./1000).^-1;
            closeAngle = par.closeAngle;
            
            % create a Fibonacci sphere to evenly distribute points
            samplingFactor = round(400/dx);      % change this if you need more points
            samplingScale = round((4*pi*r^2)/(4*pi)*0.01);
            if samplingScale == 0
                samplingScale = 1;
            end
            samples = 1 * samplingFactor * samplingScale;
            rnd = 1;
            offset = 2./samples;
            increment = pi * (3 - 5^0.5);
            
            k = 1:samples;
            y = (k+0.5).* offset - 1;
            
            rs = zeros(size(k));
            rs(1:(round(samples/2))) = (1 - y(1:(round(samples/2))).^2).^0.5;
            rs(round(samples/2)+1:samples) = rs((floor(samples/2)):-1:1);          % reduced the computational cost by 
            phi = rem((k + rnd), samples) .* increment;

            x = cos(phi) .* rs(k);
            z = sin(phi) .* rs(k);
            
            z = real(z);
            indKept = z < sin(closeAngle*pi/180);
            model.y = real(r*y(indKept))';
            model.x = real(r*x(indKept))';
            model.z = -r*z(indKept)';
            p.zOffset = median(model.z);
            model.z = model.z - p.zOffset;
            model.n = ones(size(model.x));
        end
    end
end