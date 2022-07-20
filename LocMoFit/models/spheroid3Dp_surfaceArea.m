classdef spheroid3Dp_surfaceArea<spheroidCap3Dp_surfaceArea
    % Describing endocytic coat proteins as molecules covering a part of
    % sphere with an angle indicating the closed part.
    %
    % Last update:
    %   18.07.2021
    methods
        function obj = spheroid3Dp_surfaceArea(varargin)
            obj@spheroidCap3Dp_surfaceArea(varargin{:});
            % Define the default argument values here in the constructor.
            obj.name = {'surfaceArea', 'flattening'};
            obj.fix = [0 0] ;
            obj.value = [5 0];
            obj.lb = [-inf -inf];
            obj.ub = [inf inf];
            obj.min = [0 0];
            obj.max = [50 0.7];
            
            % Define other properties here:
            obj.modelType = 'discretized';
            obj.modelTypeOption = {'discretized', 'continuous'};
            obj.dimension = 3;
            obj.listed = true;
            
            % Specific to a parametric model
        end

        function [model, p] = definedModel(obj, u, v, par, dx)
            [model, p] = definedModel@spheroidCap3Dp_surfaceArea(obj, u, v, par, dx);
            % exchange y and z to match the convention that xy is
            % rotationally symmetric.
            y = model.y;
            z = model.z;
            model.y = z;
            model.z = y;
        end
        
        function par = convertPar(obj, par)
            % spheroid is a special form of spheroid cap when closeAngle =
            % 180
            par.closeAngle = 180;
            par = convertPar@spheroidCap3Dp_surfaceArea(obj,par);
        end

        function derivedPars = getDerivedPars(obj, par)
            par.closeAngle = 180;
            derivedPars = getDerivedPars@spheroidCap3Dp_surfaceArea(obj,par);
        end
    end
end