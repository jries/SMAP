classdef locsModel<functionModel
    properties
        pseudoModel
        linkedLocsImgModel
    end
    methods
        function obj = locsModel(locs, varargin)
            % The constructor of the functional model object. This function
            % fetches the default values from the geometric model.
            p = inputParser;
            p.addParameter('layer', 1);
            parse(p,varargin{:});
            layer = p.Results.layer;
            if isfield(locs, 'z')
                obj.dimension = 3;
            else
                obj.dimension = 2;
            end
%             obj.dimension = length(size(img));
            obj.pixelSize = 1;
            obj.modelType = 'discrete';
            % subsetting the locs according to the specified layer
            fn = fieldnames(locs);
            for k=1:length(fn)
                locs.(fn{k}) = locs.(fn{k})(locs.layer == layer);
            end
            obj.pseudoModel = locs;
            addlistener(obj, 'mParsArgModified', @mParsArgModified_callback);
        end
        function set.pseudoModel(obj, val)
            obj.pseudoModel = val;
            if ~isempty(obj.linkedLocsImgModel)
                obj.linkedLocsImgModel.cleanImg();
            end
        end
        function [ref,internalPar] = getPoint(obj, mPar, varargin)
            if ~isempty(obj.ParentObject)&&length(obj.pseudoModel)>1
                fitter = obj.ParentObject;
                version = fitter.modelVerCascade(fitter.currentCascadeStep);
            else
                version = 1;
            end
            ref.x = obj.pseudoModel(version).xnm;
            ref.y = obj.pseudoModel(version).ynm;
            % Only load znm for 3D data
            if isfield(obj.pseudoModel,'znm')
                ref.z = obj.pseudoModel(version).znm;
            end
            if isfield(obj.pseudoModel, 'n')
                ref.n = obj.pseudoModel(version).n;
            else
                ref.n = ones(size(ref.x));
            end
            internalPar = [];
        end
    end
end