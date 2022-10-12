classdef Loader4PiMat<WorkflowModules.Loaders.TifLoader
    properties
        guidefh
    end
    methods
        function obj=Loader4PiMat(varargin)
            obj@WorkflowModules.Loaders.TifLoader(varargin{:})
            obj.loaders={'mat',@imageloader_mat};
            obj.fileext='*.mat';
        end   
        function initGui(obj)
            initGui@WorkflowModules.Loaders.TifLoader(obj);
            obj.guihandles.ismultifile.Visible='off';
            obj.guihandles.onlineanalysis.Visible='off';
            obj.guihandles.onlineanalysiswaittime.Visible='off';
        end
    end
end