classdef continousLinearModel_midPoint<cspline3D_midPoint
    % This class is just a container. It is replaced by :class:`cspline3D_midPoint`
    methods
        function obj = continousLinearModel_midPoint(varargin)
            obj@cspline3D_midPoint(varargin{:}); 
            obj.listed = false;
        end
    end
end