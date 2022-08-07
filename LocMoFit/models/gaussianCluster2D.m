classdef gaussianCluster2D<geometricModel
    % :class:`sphericalCap3Dp_surfaceArea` describes a geometry of spherical cap in 3D. It describe the same geometry with the same parameterization as :class:`sphericalCap3D_surfaceArea<models.sphericalCap3D_surfaceArea>` but in a parametric form.
	%
	% Geometric parameters:
    %   * `x0`: (nm) the x position of the cluster.
    %   * `y0`: (nm) the y position of the cluster.
	%
	% .. important::
	%
	%    The parameter sigma of the gaussian cluster is determined by the extrinsic parameter `variation`.
	%
    % Relavent biological structure:
    %   * a protein cluster on the plasma membrane
	%
    % Preview:
	% 	.. image:: ./images/models/gaussianCluster2D.PNG
    %       :width: 400
    methods
        function obj = gaussianCluster2D(varargin)
            obj@geometricModel(varargin{:});
            obj.name = {'x0', 'y0'};
            obj.fix = [0 0] ;
            obj.value = [0 0];
            obj.lb = [-inf -inf];
            obj.ub = [inf inf];
			obj.min = [-50 -50];
            obj.max = [50 50];
			
            obj.modelType = 'discretized';
            obj.modelTypeOption = {'discretized'};
            obj.dimension = 2;
            obj.listed = true;
        end
        
        function [model, p]= reference(obj, par, dx)
            model.x = par.x0;
            model.y = par.y0;
            model.n = 1;
            p = [];
        end
    end
end
