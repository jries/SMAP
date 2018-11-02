classdef LayerInterface< interfaces.DialogProcessor
%     classdef LayerInterface< interfaces.GuiModuleInterface & interfaces.LocDataInterface
    %provides functionality for layers
    properties
        layer=0  %layer number
    end
    methods
        function obj=LayerInterface(varargin)    
%             obj@interfaces.GuiModuleInterface(varargin{:});
            obj@interfaces.DialogProcessor(varargin{:});
        end
        function fpref=layerprefix(obj) 
            %prefix for current layer. Use for obj.getPar
            fpref=['layer' num2str(obj.layer) '_'];
        end
%         function setLayerParameter
    end
end