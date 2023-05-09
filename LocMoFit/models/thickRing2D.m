classdef thickRing2D<geometricModel
    % :class:`thickRing2D` describes the side-view projection of a thick ring.
    %
    % Geometric parameters:
    %   * `innerRadius`: (nm) the inner radius of the ring.
    %   * `outerRadius`: (nm) the outer radius of the ring.
	%   * `thickness`: (nm) the thickness of the ring.
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
    %   .. image:: ./images/models/thickRing2D.PNG
    %       :width: 400
    %   Scale bar: 50 nm.
    methods
        function obj = thickRing2D(varargin)
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
            obj.listed = true;
        end
        function [model, p]= reference(obj,par, dx)
		% For details, see :meth:`reference`.
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