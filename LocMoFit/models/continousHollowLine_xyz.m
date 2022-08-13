classdef continousHollowLine_xyz<csplineTube3D_xyz
% This class is just a container. It is replaced by :class:`csplineTube3D_xyz`
    methods
        function obj = continousHollowLine_xyz(varargin)
            obj@csplineTube3D_xyz(varargin{:}); 
            obj.listed = false;
        end
    end
end