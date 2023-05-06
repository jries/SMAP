classdef dualEllipse3D_avgR<dualEllipse3D_avgR_discrete
    % This class is just a container. It is replaced by :class:`dualEllipse3D_avgR_discrete`
    methods
        function obj = dualEllipse3D_avgR(varargin)
            obj@dualEllipse3D_avgR_discrete(varargin{:}); 
            obj.listed = false;
        end
    end
end
