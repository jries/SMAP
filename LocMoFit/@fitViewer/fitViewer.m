classdef fitViewer<matlab.mixin.Copyable
    % :class:`fitViewer` allows the user to interactively visualize the fit result by LocMoFit. This class creates an interactive viewer for inspeting the result in 3D.
	%
	% Usage:
	%
	% Attributes:
	%	linkedLocMoFit (LocMoFit object)
    %	localizations (cell)
    %	modelPoints (cell)
    %
	% Last update:
    %   26.12.2022
    %
    % See also:
    %   :class:`LocMoFit<@LocMoFit.LocMoFit>`

    properties
        linkedLocMoFit
        localizations
        modelPoints
    end
    methods
        function obj = fitViewer(varargin)
        end

        function prepareData(obj)

        end

        function imgRendering(obj)

        end

        function createGUI(obj)
            % create the GUI

            obj.GUILayout

            obj.freeRotView
            
            obj.fixView
        end

        function freeRotView(obj)
            % create a view in which the data can be rotated
        end
        
        function fixView(obj)
            % create a view fixed at a certain orientation
        end

        function GUILayout(obj)
        end

        function interaction(obj)

        end
    end
end