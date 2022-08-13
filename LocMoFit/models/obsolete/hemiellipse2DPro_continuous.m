classdef hemiellipse2DPro_continuous<hemispheroid2D
    % This class is just a container. It is replaced by :class:`hemispheroid2D`
    %
    % Last update:
    %   21.07.2021
    methods
        function obj = hemiellipse2DPro_continuous(varargin)
            obj@hemispheroid2D(varargin{:}); 
            obj.listed = false;
        end
    end
end