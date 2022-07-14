classdef geometricModel<matlab.mixin.Copyable
    % :class:`geometricModel` is the superclass of any geometric model. It
    % contains methods for building own geometric models.
    %
    % Last update:
    %   14.10.2021

    properties (SetObservable)
        ParentObject = [];
        name                 % names of model parameters
        fix                  % fixing model parameters
        value                % values of model parameters
        lb                   % relative lower bounds of model parameters
        ub                   % relative upper bounds of model parameters
        min                  % min values of model parameters
        max                  % max values of model parameters
        internalSettings     % parameters that do not suit fitting
    end
    properties
        modelType            % selected model type
        modelTypeOption      % possible model types of a specific geometric model
        dimension            % a scalar indicating the dimensionality of the model. Either 2 or 3.
        listed = false       % whether this model will be listed in the GUI.
    end
    properties (Hidden)
        parsArgName = {'name', 'fix', 'value', 'lb', 'ub', 'min', 'max'}';
    end
    methods
        function obj = geometricModel(varargin)
            % Usage:
            %   obj = geometricModel(varargin)
            %
            % Args:
            %   Name-value pairs:
            %   * 'Parent': parental obj.
            %
            % Returns:
            %   obj: an :class:`geometricModel` object.
                
            if length(varargin)>0
                p = inputParser;
                p.addParameter('Parent', [])
                parse(p,varargin{:});
                obj.ParentObject = p.Results.Parent;
            end
            if isempty(obj.listed)
                obj.listed = false;
            end
        end
        function setInternalSettings(obj,settingName,value)
            if isfield(obj.internalSettings, settingName)
                obj.internalSettings.(settingName) = value;
            else
                warning([settingName ' is not a internal property of ' class(obj) '.']);
            end
        end
        function value = getInternalSettings(obj,settingName)
            if isfield(obj.internalSettings, settingName)
                value = obj.internalSettings.(settingName);
            else
                warning([settingName ' is not a internal property of ' class(obj) '.']);
            end
        end
        function value = internalSettingsModified_callback(obj,event)
            value = [];
            disp('No internalSettingsModified_callback defined.');
        end
        %% Model parameter-related methods
        function addMPar(obj,newMPar)
            parsArgName = obj.parsArgName;
            n = length(newMPar.name);
            for k = 1:length(parsArgName)
                obj.(parsArgName{k})(end+1:end+n) = newMPar.(parsArgName{k});
            end
        end
        function rmMPar(obj,parsName)
            parsArgName = obj.parsArgName;
            n = length(parsName);
            for l = 1:n
                idx = strcmp(obj.name, parsName{l});
                for k = 1:length(parsArgName)
                    obj.(parsArgName{k})(idx) = [];
                end
            end
        end
        
        function items = getThings2Plot(obj,mPar)
        % The user can define what should be also displayed in the plots.
        % Args:
        % 	obj: a :class:`geometricModel` object.
        % Returns:
        %   items: things to be plotted.
        %
            items = [];
        end
        
        function derivedPars = getDerivedPars(obj, varargin)
        % Get derived parameters of the current model.
        % Args:
        % 	obj: an :class:`functionModel` object.
        % Returns:
        %   derivedPars: derived parameters.
        %
            derivedPars = [];
        end
        
        %% Continuous-to-discrete conversion
%         function modelRef = con2dis(obj, modelRef)
%             lDuplicate = any(modelRef.n>1);
%             if lDuplicate
%                 find(modelRef.n>1);
%                 
%             else
%                 modelRef = rmfield(modelRef,'n');
%             end
%         end
        
        %% Response when set
        function set.internalSettings(obj, value)
            obj.internalSettings = [];
            if isstruct(value)
                fn = fieldnames(value);
                for k = 1:length(fn)
                    obj.internalSettings.(fn{k}) = value.(fn{k});
                end
                respond2InternalSettingsChange(obj);
            else
                obj.internalSettings = value;
            end
        end
        
    end
    methods(Access = protected)
        function respond2InternalSettingsChange(obj)
        end
    end
end
