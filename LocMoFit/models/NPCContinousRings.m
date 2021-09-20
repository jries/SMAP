classdef NPCContinousRings
    properties
        name = {'ringDistance', 'radius'};
        fix = [0 1] ;
        value = [0 53.7];
        min = [0 0];
        max = [90 70];
        lb = [0 0];
        ub = [90 0];
        modelType = 'continuous';
        dimension = 3;
    end
    methods
        function obj = NPCContinousRings
        end
    end
    methods (Static)
        function [model, p]= reference(par, dx)
            r = par.radius;
            d = par.ringDistance;

            numphi=round(2*pi*r/dx);
            theta = linspace(0,2*pi,numphi);
            theta(end)=[];
            model = [];
            x = r.*cos(theta); y = r.*sin(theta);
            z = repelem(d/2, length(theta));

            model.x = repmat(x,1,2); model.y = repmat(y,1,2);
            model.z = [z -z];
            model.n = ones(size(model.x));
            p = [];
        end
    end
end