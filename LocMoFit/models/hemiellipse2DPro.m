classdef hemiellipse2DPro
    properties
        name = {'xcenter', 'ycenter', 'a', 'b'};
        fix = [1 1 1 1] ;
        value = [0 0 85*0.9 40];
        lb = [-60 -230 -10 -10];
        ub = [60 30 10 10];
        modelType = 'continuous'
        dimension = 2;
        min = [-150 -150 10 10]
        max = [150 150 40 200]
    end
    methods
        function obj = hemiellipse2DPro
        end
    end
    methods (Static)
        function [model, p]= reference(par, dx)
            c = [par.xcenter par.ycenter par.a par.b];
            roiSize = 300;
            num = round(roiSize/dx)*6;
            [model.x, model.y] = meshgrid(linspace(-roiSize/2,roiSize/2,num), linspace(-roiSize/2,roiSize/2,num));
            model.x = model.x(:);
            model.y = model.y(:);
            model.n = hemiellipse(c,model.x, model.y);
            p = [];
        end
    end
end

function v = hemiellipse(c,x,y)
% c = [xcenter, ycenter, a, b]
% x and y are the coordinates of the grid evenly distributed in ROI
    r = x - c(1); % distance to xcenter
    l = y - c(2); % distance to ycenter
    
    % in y-dim
    rCircles = ellipse(c(3),c(4),l);
    
    l = l >= -c(4)& l <= 0; % set l to zero if the coordinate is out of range
    
    % in x-dim
    v = disk(rCircles, r);
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

function v = ellipse(a,b,d)
    % a: half horizental axis
    % b: half vertical axis
    % d: distance to the center
    v = zeros(size(d));
    idx = abs(d)<abs(b);
    if length(d)>1
        idxD = idx;
    else
        idxD = 1;
    end
    if any(idx)
        v(idxD) = a.*cos(asin(d(idxD)./b));
    end
end
