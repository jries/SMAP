classdef spheroid3Dp_surfaceArea<spheroidCap3Dp_surfaceArea
    % :class:`spheroid3Dp_surfaceArea` describes the geometry of spheroid in 3D. It is parametric. Spheroid is a sphere flattened at the poles. 
	%
	% .. Important::
	%
	%   Here the flattening is applied along the z-axis, not the y-axis.
	%
	% Geometric parameters:
    %   * `surfaceArea`: (10\ :sup:`4` nm\ :sup:`2`) the surface area of the spherical cap.
	%	* `flattening`: (no unit) or `f`, is defined as `1-c/a`, where a and c are the two distinct axis lengths. c lines on the y-axis.
	%
    % Relavent biological structure:
    %   * deformed vesicle
	%
	% See also:
    %   :class:`spheroidCap3Dp_surfaceArea<models.spheroidCap3Dp_surfaceArea>`
	%
    % Preview:
	% 	.. image:: ./images/models/spheroid3Dp_surfaceArea.PNG
    %       :width: 400
    %   Scale bar: 50 nm.
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