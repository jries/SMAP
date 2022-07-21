classdef thickRing2DPro_continuous<thickRing2D
    % This class is just a container. It is replaced by :class:`thickRing2D`
    %
    % Last update:
    %   21.07.2021
    methods
        function obj = thickRing2DPro_continuous(varargin)
            obj@thickRing2D(varargin{:}); 
            obj.listed = false;
        end
    end
end