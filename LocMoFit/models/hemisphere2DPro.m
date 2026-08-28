classdef hemisphere2DPro<geometricModel
    % :class:hemisphere2DPro describes the side-view projection of a hemisphere.
    %
    % Geometric parameters:
    %   * `xcenter`: (nm) [obsolete] please set it to zero.
    %   * `ycenter`: (nm) [obsolete] please set it to zero.
    %   * `innerRadius`: (nm) radius of the inner ring.
    %   * `outerRadius`: (nm) radius of the outer ring.
    %
    % Relevant biological structure:
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
        function obj = hemisphere2DPro(varargin)
            obj@geometricModel(varargin{:});
            obj.name = {'xcenter', 'ycenter', 'innerRadius', 'outerRadius'};
            obj.fix = [1 1 1 1] ;
            obj.value = [0 0 20 40];
            obj.lb = [-60 -230 -10 -10];
            obj.ub = [60 30 10 10];
            obj.min = [-60 -230 -10 -10];
            obj.max = [60 30 10 10];
            
            % Define other properties here:
            obj.modelType = 'continuous';
            obj.modelTypeOption = {'discretized' 'continuous'};
            obj.dimension = 2;
            obj.listed = true;
        end
        
        function [model, p]= reference(obj, par, dx)
		% For details, see :meth:`reference`.
            c = [par.xcenter par.ycenter par.innerRadius par.outerRadius];
            fitter = obj.ParentObject.ParentObject;
            
            bond = fitter.roiSize/2 + fitter.imgExtension;
            
            if isempty(obj.ParentObject.locsPrecFactor)
                locsPrecFactor = 1;
            else
                locsPrecFactor = obj.ParentObject.locsPrecFactor;
            end
            
            num = round(bond*2/(dx*locsPrecFactor));
            [model.x, model.y] = meshgrid(linspace(-fitter.roiSize/2,fitter.roiSize/2,num), linspace(-fitter.roiSize/2,fitter.roiSize/2,num));
            model.x = model.x(:);
            model.y = model.y(:);
            model.n = cap(c,model.x, model.y);
            p = [];
        end
        
        function derivedPars = getDerivedPars(obj, pars)
            derivedPars = [];
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
