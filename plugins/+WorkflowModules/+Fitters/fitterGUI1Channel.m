classdef fitterGUI1Channel<WorkflowModules.Fitters.fitterGUI
    %     Intermediate GUI to select a fitting plugin.
    methods
        function obj=fitterGUI1Channel(varargin)
            obj@WorkflowModules.Fitters.fitterGUI(varargin{:})
            obj.inputChannels=1; %1: image. 2: background image
        end
    end
end