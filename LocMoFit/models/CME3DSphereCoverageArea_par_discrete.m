classdef CME3DSphereCoverageArea_par_discrete<sphericalCap3Dp_surfaceArea
    % This class is just a container. It is replaced by :class:`sphericalCap3D_surfaceArea_par`
    methods
        function obj = CME3DSphereCoverageArea_par_discrete(varargin)
            obj@sphericalCap3Dp_surfaceArea(varargin{:}); 
            obj.listed = false;
        end
    end
end