classdef CME3DSphereCoverageArea_discrete<sphericalCap3D_surfaceArea
    % This class is just a container. It is replaced by :class:`sphericalCap3D_surfaceArea`
    methods
        function obj = CME3DSphereCoverageArea_discrete(varargin)
            obj@sphericalCap3D_surfaceArea(varargin{:}); 
            obj.listed = false;
        end
    end
end