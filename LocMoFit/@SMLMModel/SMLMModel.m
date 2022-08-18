classdef SMLMModel<matlab.mixin.Copyable
    % :class:`SMLMModel` is a super-class for defining how to deal with a geometric
    % model.
    %
    % Todo:
    %   *The modelType needs to be extended to discretized and intensity.
    properties (SetObservable)
        ParentObject = [];      % Parental SMLMModelFit object.
        ID                      % The model's ID in the parental SMLMModelFit object.
        img                     % Model image.
        parVal                  % [obsolete].
        mPars                   % Model parameters.
        modelObj                % Source geometric model object.
        modelFun                % The function for creating coordinates based on the geometric model.
        sourcePath              % The path of the m file of the geometric model.
        dimension               % Dimension of the geometric model.
        modelType               % Type of the model, either discrete, discretized, continuous, intensity, or image.
        weight = 1;
        fixSigma = false;       % Fix the sigma to a specific value.
        displayLut = 'red hot'; % The lookup table for the model.
        layer = 1;              % The layer that this model is fitted to.
    end
    
    methods
        function obj = SMLMModel()
        end
        function addParent(obj, parent)
            % Add the parental SMLMModelFit object.
            obj.ParentObject = parent;
        end
        intensityVal = modelHandler(obj, locs, mPars);
        function mPars = exportMPars(obj)
            % Export model parameters and their default values.
            for k = 1:length(obj.mPars.name)
                mPars.(obj.mPars.name{k})=obj.mPars.value(k);
            end
        end
     function set.modelObj(obj,val)
           obj.modelObj = val; % updates the property
           respond2ModelObjChange(obj) % calls the protected method
        end 
     end 
     methods(Access = protected)
        function respond2ModelObjChange(obj, value) % no function definition here
        end
     end
    events
           mParsArgModified      % When the list of model parameters being changed.
    end
end
function likelihood = gaussDist(refLocs,y,x,z,ysigmao,xsigmao,zsigmao, varargin)
    p = inputParser;
    addParameter(p, 'sigFactor',1);
    parse(p, varargin{:})
    % For one dimension: exp(-((x-mu/sigma)^2)/2)/(2*pi)^(1/2)*sigma
    f = p.Results.sigFactor;          % scale factor of the 3 sigma
    xsigma = xsigmao.*f;ysigma = ysigmao.*f;zsigma = zsigmao.*f;
    
    % the template bellow is for multiple sets of parameters at the same
    % time
    templateAll = zeros([1 size(refLocs.x,2) size(refLocs.x,1)]);
    template.x = templateAll;
    template.y = templateAll;
    template.z = templateAll;
    template.x(1,:,:)=refLocs.x';
    template.y(1,:,:)=refLocs.y';
    template.z(1,:,:)=refLocs.z';
    refLocs = template;
    termX = ((x - refLocs.x)./xsigma).^2;
    termY = ((y - refLocs.y)./ysigma).^2;
    termZ = ((z - refLocs.z)./zsigma).^2;
    svol = ((2*pi)^(-3/2))./(xsigma.*ysigma.*zsigma);
    likelihood = sum(svol.*exp(-(termX+termY+termZ)/2),3)';
    likelihood = likelihood';
end