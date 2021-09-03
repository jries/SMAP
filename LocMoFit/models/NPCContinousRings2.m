classdef NPCContinousRings2<geometricModel
    methods
        function obj = NPCContinousRings2(varargin)
            obj@geometricModel(varargin{:});
            obj.name = {'ringDistance', 'radius'};
            obj.fix = [0 1] ;
            obj.value = [0 53.7];
            obj.min = [0 0];
            obj.max = [90 70];
            obj.lb = [0 0];
            obj.ub = [90 0];
            obj.modelType = 'continuous';
            obj.dimension = 3;
        end
        
        function [model, p]= reference(obj, par, dx)
            r = par.radius;
            d = par.ringDistance;

            numphi=round(2*pi*r/dx);
            theta = linspace(0,2*pi,numphi);
            theta(end)=[];
            model = [];
            x = r.*cos(theta); y = r.*sin(theta);
            z = repelem(d/2, length(theta));

            model.x = repmat(x,1,2)'; model.y = repmat(y,1,2)';
            model.z = [z -z]';
            model.n = ones(size(model.x));
            p = [];
        end
    end
end