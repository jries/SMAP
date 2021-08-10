classdef simplex2D<geometricModel
    % simplex2D generates the vertex of a triangle, the simplex in 2D.
    
    methods
        function obj = simplex2D(varargin)
            obj@geometricModel(varargin{:});
            obj.name = {'edgeLength'};
            obj.fix = 1;
            obj.value = 15;
            obj.lb = 0;
            obj.ub = 0;
            obj.modelType = 'discrete';
            obj.dimension = 2;
            obj.min = 15;
            obj.max = 15;
            obj.internalSettings.replicate = 5;
        end
        
        function [model, p]= reference(obj, par, dx)
            edgeLength = par.edgeLength;
            model.x = repmat([-edgeLength/2 0 edgeLength/2]',obj.internalSettings.replicate,1);
            model.y = repmat([-edgeLength/2 edgeLength/2 -edgeLength/2]',obj.internalSettings.replicate,1);
            model.n = ones(obj.internalSettings.replicate*3,1);
            p = [];
        end
    end
end
