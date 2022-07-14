classdef CME3DSphereCoverageArea_par_discrete<sphericalCap3D_surfaceArea_par
    % This class is just a container. It is replaced by :class:`sphericalCap3D_surfaceArea_par`
    methods
        function obj = CME3DSphereCoverageArea_par_discrete(varargin)
            obj@sphericalCap3D_surfaceArea_par(varargin{:}); 
            obj.listed = false;
        end
    end
end