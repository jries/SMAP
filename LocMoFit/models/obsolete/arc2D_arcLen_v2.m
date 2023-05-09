classdef arc2D_arcLen_v2<arc2D_arcLen
    % This class is just a container. It is replaced by :class:`arc2D_arcLen`
    %
    % Last update:
    %   21.07.2021
    methods
        function obj = arc2D_arcLen_v2(varargin)
            obj@arc2D_arcLen(varargin{:}); 
            obj.listed = false;
        end
    end
end