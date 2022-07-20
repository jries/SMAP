classdef NPCPointModel_flexible2<dualRing3D_discrete
    % This class is just a container. It is replaced by :class:`dualRing3D_discrete`
    methods
        function obj = NPCPointModel_flexible2(varargin)
            obj@dualRing3D_discrete(varargin{:}); 
            obj.listed = false;
        end
    end
end
