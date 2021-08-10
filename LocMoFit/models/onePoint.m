classdef onePoint<geometricModel
    methods
        function obj = onePoint(varargin)
            obj@geometricModel(varargin{:});
            obj.name = {};
            obj.fix = [] ;
            obj.value = [];
            obj.lb = [];
            obj.ub = [];
            obj.modelType = 'discrete';
            obj.dimension = 3;
            obj.min = [];
            obj.max = [];
        end
        
        function [model, p]= reference(obj, par, dx)
            model.x = 0;
            model.y = 0;
            model.z = 0;
            model.n = 1;
            p = [];
        end
    end
end
