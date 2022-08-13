classdef dualRingModel<dualRing3D_discrete
    % This class is just a container. It is replaced by :class:`dualRing3D`
    methods
        function obj = dualRingModel(varargin)
            obj@dualRing3D_discrete(varargin{:}); 
            obj.listed = false;
        end
    end
end
