classdef LocDataInterface<handle;
    %provides access to localization Data
    properties
        locData %interfaces.LocalizationData object, here localizations are saved   
    end

    methods
        function obj=LocDataInterface(varargin)
            if nargin>0
                obj.attachLocData(varargin{1}); 
            else
                obj.locData=interfaces.LocalizationData();
            end
        end
        function attachLocData(obj,locData)
            %sets obj.locData to handle object locData, to share it among
            %all modules
            obj.locData=locData;
        end

    end
end