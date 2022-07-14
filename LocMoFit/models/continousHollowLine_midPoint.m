classdef continousHollowLine_midPoint<csplineTube3D_midPoint
% This class is just a container. It is replaced by :class:`csplineTube3D_midPoint`
    methods
        function obj = continousHollowLine_midPoint(varargin)
            obj@csplineTube3D_midPoint(varargin{:}); 
            obj.listed = false;
        end
    end
end