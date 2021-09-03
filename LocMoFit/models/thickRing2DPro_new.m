classdef thickRing2DPro_new<geometricModel
    methods
        function obj = thickRing2DPro_new(varargin)
            obj@geometricModel(varargin{:});
            obj.name = {'xcenter', 'ycenter', 'innerRadius', 'outerRadius', 'thickness'};
            obj.fix = [1 1 1 1 1] ;
            obj.value = [0 0 40 70 60];
            obj.lb = [-60 30 -10 -10 -5];
            obj.ub = [60 230 10 10 5];
            obj.min = [-60 30 20 80 40];
            obj.max = [60 230 80 100 70];
            obj.modelType = 'continuous';
            obj.dimension = 2;
        end
    end
    methods
        function [model, p]= reference(obj,par, dx)
            c = [par.xcenter par.ycenter par.innerRadius par.outerRadius par.thickness];
            roiSize = 300;
            num = round(roiSize/dx)*6;
            [model.x, model.y] = meshgrid(linspace(-roiSize/2,roiSize/2,num), linspace(-roiSize/2,roiSize/2,num));
            model.x = model.x(:);
            model.y = model.y(:);
            model.n = thickRing(c,model.x, model.y);
            p = [];
        end
    end
end

function v = thickRing(c,x,y)
% c = [xcenter, ycenter, innerDia, outerDia, thickness]
    r = x - c(1); 
    l = y - c(2);
    
    l = l >= -c(5)/2 & l <= c(5)/2;

    v = disk(c(4), r) - disk(c(3), r);
    v = v .* l;
end

function v = disk(r,d)
    % r: radius
    % d: distance to the center
    v = zeros(size(d));
    idx = abs(d)<abs(r);
    if length(r)>1
        idxD = idx;
    else
        idxD = 1;
    end
    if any(idx)
        sqv = r(idxD).^2-d(idx).^2; % based on pythagorean theorem
        v(idx) = 2.*sqrt(sqv);
    end
end