classdef triangle<geometricModel
    % simplex2D generates the vertex of a triangle, the simplex in 2D.
    
    methods
        function obj = triangle(varargin)
            obj@geometricModel(varargin{:});
            obj.name = {'edgeLength'};
            obj.fix = 1;
            obj.value = 15;
            obj.lb = 0;
            obj.ub = 0;
            obj.modelType = 'discrete';
            obj.modelTypeOption = {'discrete'};
            obj.dimension = 3;
            obj.min = 15;
            obj.max = 15;
        end
        
        function [model, p]= reference(obj, par, dx)
            edgeLength = par.edgeLength;
            model.x = [-edgeLength/2 0 edgeLength/2]';
            model.y = [-edgeLength/2 edgeLength/2 -edgeLength/2]';
            model.z = zeros(3,1);
            model.n = ones(3,1);
            p = [];
        end
    end
end
