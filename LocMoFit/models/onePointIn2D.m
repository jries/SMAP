classdef onePointIn2D<geometricModel
    methods
        function obj = onePointIn2D(varargin)
            obj@geometricModel(varargin{:});
            obj.name = {};
            obj.fix = [] ;
            obj.value = [];
            obj.lb = [];
            obj.ub = [];
            obj.modelType = 'discrete';
            obj.modelTypeOption = {'discrete'};
            obj.dimension = 2;
            obj.min = [];
            obj.max = [];
        end
        
        function [model, p]= reference(obj, par, dx)
            model.x = 0;
            model.y = 0;
            model.n = 1;
            p = [];
        end
    end
end
