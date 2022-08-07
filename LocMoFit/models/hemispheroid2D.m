classdef hemispheroid2D<geometricModel
    % :class:`hemispheroid2D` describes the side-view projection of a hemispheroid.
    %
    % Geometric parameters:
    %   * `a`: (nm) the axis along the y-axis.
    %   * `b`: (nm) the axis along the x-axis.
    %   * `xcenter`: (nm) [obsolete] please set it to zero.
    %   * `ycenter`: (nm) [obsolete] please set it to zero.
    %
    % Relavent biological structure:
    %   * actin network at the endocytic site
    %
    % Preview:
	% 	.. note::
	% 		It will be available soon.
	%
	% ..
    %   .. image:: ./images/models/hemispheroid2D.PNG
    %       :width: 400
    %   Scale bar: 50 nm.
    methods
        function obj = hemispheroid2D(varargin)
            obj@geometricModel(varargin{:});
            obj.name = {'xcenter', 'ycenter', 'a', 'b'};
            obj.fix = [1 1 1 1] ;
            obj.value = [0 0 85*0.9 40];
            obj.lb = [-60 -230 -10 -10];
            obj.ub = [60 30 10 10];
            obj.min = [-150 -150 10 10];
            obj.max = [150 150 40 200];
            
            % Define other properties here:
            obj.modelType = 'discretized';
            obj.modelTypeOption = {'discretized' 'continuous'};
            obj.dimension = 2;
            obj.listed = true;
        end
        
        function [model, p]= reference(obj, par, dx)
		% For details, see :meth:`reference`.
            c = [par.xcenter par.ycenter par.a par.b];
            fitter = obj.ParentObject.ParentObject;
            
            bond = fitter.roiSize/2 + fitter.imgExtension;
            
            if isempty(obj.ParentObject.locsPrecFactor)
                locsPrecFactor = 1;
            else
                locsPrecFactor = obj.ParentObject.locsPrecFactor;
            end
            
            num = round(bond*2/(dx*locsPrecFactor));
            [model.x, model.y] = meshgrid(linspace(-bond, bond, num), linspace(-bond, bond, num));
            model.x = model.x(:);
            model.y = model.y(:);
            model.n = hemiellipse(c,model.x, model.y);
            p = [];
        end
        
        function derivedPars = getDerivedPars(obj, pars)
            derivedPars = [];
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
