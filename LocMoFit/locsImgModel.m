classdef locsImgModel<imageModel
    properties
        linkedLocsModel
        tempImg
        sigmaFactor = [1 0];
    end
    methods
        
        function obj = locsImgModel(locs,varargin)
            % The constructor of the functional model object. This function
            % fetches the default values from the geometric model.
            % Fetches the default values
            % parse varargin
            obj@imageModel();
            
            p = inputParser;
            p.addParameter('pixelSize', 2);
            p.addParameter('layer', 1);
            parse(p,varargin{:});
            results = p.Results;
            pixelSize = results.pixelSize;
            layer = results.layer;
            obj.pixelSize = pixelSize;
%             obj.img = img;
%             obj.dimension = length(size(img));

            % subsetting the locs according to the specified layer
            fn = fieldnames(locs);
            for k=1:length(fn)
                locs.(fn{k}) = locs.(fn{k})(locs.layer == layer);
            end

            % Synchornized the linkedLocsModel
            obj.linkedLocsModel = locsModel(locs, 'layer', layer);
            obj.linkedLocsModel.pixelSize = obj.pixelSize;
%             obj.linkedLocsModel.eps = obj.eps;
            obj.linkedLocsModel.sigma = obj.gaussSigma;
            obj.linkedLocsModel.layer = obj.layer;
            obj.linkedLocsModel.linkedLocsImgModel = obj;
            
            addlistener(obj, 'mParsArgModified', @mParsArgModified_callback);
            addlistener(obj,'pixelSize','PostSet',@obj.handlePropEvents);
%             addlistener(obj,'eps','PostSet',@obj.handlePropEvents);
            addlistener(obj,'gaussSigma','PostSet',@obj.handlePropEvents);
            addlistener(obj,'layer','PostSet',@obj.handlePropEvents);
        end
        
        %% image related functions
        function loadImg(obj,varargin)
            % Load an image from the linked locsModel.
            
            % Here set the roiSize according to the parental SMLMModelFit
            % obj.
            lRoiSize = strcmp(varargin, 'roiSize');
            if sum(lRoiSize)==0
                varargin = [varargin {'roiSize', obj.ParentObject.roiSize}];
            end

            lUseLocP = strcmp(varargin, 'useLocprecnm');
            if sum(lUseLocP)==0
                varargin = [varargin {'useLocprecnm', ~obj.fixSigma}];
            end
            
            obj.linkedLocsModel.sigmaFactor = obj.sigmaFactor;
            obj.img = obj.linkedLocsModel.getImage([],varargin{:});
        end
        
        function checkImg(obj)
            if length(obj.sigmaFactor) == 2 && isfield(obj.tempImg, ['p' num2str(obj.sigmaFactor(2))])
                obj.img = obj.tempImg.(['p' num2str(obj.sigmaFactor(2))]);
            else
                obj.loadImg;
                obj.saveImg;
            end
        end
        
        function saveImg(obj)
            obj.tempImg.(['p' num2str(obj.sigmaFactor(2))]) = obj.img;
        end
        
        function cleanImg(obj)
            obj.tempImg = [];
        end
        
        function set.linkedLocsModel(obj, val)
            obj.cleanImg();
            obj.linkedLocsModel = val; 
        end
        
        function handlePropEvents(obj,src,bb)
            switch src.Name
                case {'layer','pixelSize'}
                    obj.linkedLocsModel.(src.Name) = obj.(src.Name);
                case 'gaussSigma'
                    obj.linkedLocsModel.sigma = obj.(src.Name);
            end
        end
    end
end