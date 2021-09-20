classdef thickRing2DPro_continuous<geometricModel
    methods
        function obj = thickRing2DPro_continuous(varargin)
            obj@geometricModel(varargin{:});
            obj.name = {'xcenter', 'ycenter', 'innerRadius', 'outerRadius', 'thickness'};
            obj.fix = [1 1 1 1 1] ;
            obj.value = [0 0 40 70 60];
            obj.lb = [-60 30 -10 -10 -5];
            obj.ub = [60 230 10 10 5];
            obj.min = [-60 30 20 80 40];
            obj.max = [60 230 80 100 70];
            
            % Define other properties here:
            obj.modelType = 'discretized';
            obj.modelTypeOption = {'discretized' 'continuous'};
            obj.dimension = 2;
        end
        function [model, p]= reference(obj,par, dx)
            c = [par.xcenter par.ycenter par.innerRadius par.outerRadius par.thickness];
            fitter = obj.ParentObject.ParentObject;

            bond = fitter.roiSize/2 + fitter.imgExtension;
            
            if isempty(obj.ParentObject.locsPrecFactor)
                locsPrecFactor = 1;
            else
                locsPrecFactor = min(obj.ParentObject.locsPrecFactor,15);
            end
            
            num = round(bond*2/(dx*locsPrecFactor));
            [model.x, model.y] = meshgrid(linspace(-bond,bond,num), linspace(-bond,bond,num));
            model.x = model.x(:);
            model.y = model.y(:);
            model.n = thickRing(c,model.x, model.y);
            p = [];
        end
        function derivedPars = getDerivedPars(obj, pars)
            derivedPars = [];
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