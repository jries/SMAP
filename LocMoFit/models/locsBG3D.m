classdef locsBG3D<geometricModel
    
    % log
    %   - 201229: change the sign of the ring twist
    methods
        function obj = locsBG3D(varargin)
            obj@geometricModel(varargin{:});

            % Define other properties here:
            obj.modelType = 'background';
            obj.modelTypeOption = {'background','continuous'};
            obj.dimension = 3;
        end
        
        function [model, p]= reference(obj, par, dx)
            
        % Sample coordinates of the model as reference.
        % --- Syntax ---
        % [model, p]= reference(obj, par, dx)
        % --- Arguments ---
        % -- Input --
        % obj:
        % par: a structure object. Its fieldnames should be the names of
        % parameters, and their correspoinding content should be the
        % parameter values.
        % dx: sampling rate.
        % -- Output --
        % model: a structure object. Its fieldnames should be x, y, z, and
        % n, indicating the xyz position amplitude n of the sampled model
        % points.
        % p: additional information of the model.
        
        locs = obj.ParentObject.ParentObject.locs;
        lPar = obj.ParentObject.ParentObject.exportPars(1,'lPar');
        locs = obj.ParentObject.ParentObject.locsHandler(locs,lPar);
        model.x = locs.xnm-obj.ParentObject.pixelSize;
        model.y = locs.ynm-obj.ParentObject.pixelSize;
        model.z = locs.znm-obj.ParentObject.pixelSize;
        model.channel = ones(size(model.x));
        model.n = ones(size(model.x));

        %     likelihood = gaussDist(model,double(locs.ynmrot), double(locs.xnmrot), double(locs.znm), 8,8,10);
        %     likelihoodBoun = gaussDist(model,0,0,0, 16,16,20);

        p = [];
        end
        function derivedPars = getDerivedPars(obj, pars)
            derivedPars = [];
        end
    end
end
