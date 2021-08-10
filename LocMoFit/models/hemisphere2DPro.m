classdef hemisphere2DPro
    properties
        name = {'xcenter', 'ycenter', 'innerRadius', 'outerRadius'};
        fix = [1 1 1 1] ;
        value = [0 0 20 40];
        lb = [-60 -230 -10 -10];
        ub = [60 30 10 10];
        min = [-60 -230 -10 -10];
        max = [60 30 10 10];
        modelType = 'continuous'
        dimension = 2;
    end
    methods
        function obj = hemisphere2DPro
        end
    end
    methods (Static)
        function [model, p]= reference(par, dx)
            c = [par.xcenter par.ycenter par.innerRadius par.outerRadius];
            roiSize = 300;
            num = round(roiSize/dx)*6;
            [model.x, model.y] = meshgrid(linspace(-roiSize/2,roiSize/2,num), linspace(-roiSize/2,roiSize/2,num));
            model.x = model.x(:);
            model.y = model.y(:);
            model.n = cap(c,model.x, model.y);
            p = [];
        end
    end
end

function v = cap(c,x,y)
% c = [xcenter, ycenter, innerRadius, outerRadius]
% x and y are coordinates of the evenly distributed grid points
% the v means the value at y position, given a set of x-coordinate (radius)
    r = x - c(1); % distance to xcenter
    l = y - c(2); % distance to ycenter
    
    % in y-dim
    lout = disk(c(4), l);
    lin = disk(c(3), l);
    
    l = l >= -c(4)& l <= 0; % set l to zero if the coordinate is out of range
    
    % in x-dim
    v = disk(lout/2, r) - disk(lin/2, r);
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
