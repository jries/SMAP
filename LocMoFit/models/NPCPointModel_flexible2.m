classdef NPCPointModel_flexible2<dualRingModel
    % This class is just a container. It is replaced by :class:`dualRingModel`
    methods
        function obj = NPCPointModel_flexible2(varargin)
            obj@dualRingModel(varargin{:}); 
            obj.listed = false;
        end
    end
end
