classdef LocCollector<interfaces.WorkflowModule
%    collects locs
    properties
        locdata     
    end
    methods
       function obj=LocCollector(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
            obj.inputChannels=1; 
            obj.setInputChannels(1,[],'fitted localizations');
       end
        function pard=guidef(obj)
            pard.plugininfo.type='WorkflowModule'; 
            pard.plugininfo.description='';         
        end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
        end
        function prerun(obj,p)
            obj.locdata=interfaces.LocalizationData;
        end
        function output=run(obj,data,p)
            if ~isempty(data.data)    
            obj.locdata.addLocData(data.data);
            end
            output=data;
        end

    end
end

