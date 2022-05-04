classdef functionModel<SMLMModel
    % A sub-class of :class:`SMLMModel`. `functionModel` class handles any
    % geometric model in the form of a function.
    % :class:`functionModel` handles the function differently based on its
    % `modelType`. The `modelType` is per geometric model and defined in
    % `modelType` of the :class:`geometricModel`.
    % 
    % Last update:
    %   14.10.2021
    %
    % See also:
    %   :class:`SMLMModel`, :class:`LocMoFit`, :class:`geometricModel`
    
    properties
        pixelSize = 5;          % Pixel size of the model
        sigma = 15;             % Standard deviation of the gaussian kernel used for smoothing the model.
        sigmaFactor = 1;        % The scaling factor of the kernel's standard deviation.
        samplingFactor = 0.75;  % For continuous model, deciding the distance between ref points. In the unit of sigma. 0.75 means 0.75*sigma.
        sigmaZFactor = 3;
        sigmaSet                % The set of sigma.
        sigmaZSet               % The set of sigma in Z.
        extraBlurr              % This is a parameter determined by lPars.variation.
    end
    properties (Dependent)
        locsPrecFactor          % The min sqrt(locprec^2+varation^2)
    end
    methods
        intensityVal = modelHandler(obj, locs, mPars,varargin)
        function obj = functionModel(filePath)
            % The constructor of the functional model object. This function
            % fetches the default values from the geometric model.
            if exist('filePath','var')
                obj.sourcePath = filePath;
                [filePath, modelFun, ext] = fileparts(filePath);
                addpath(filePath)
                modelObj = str2func(modelFun);
                modelObj = modelObj('Parent', obj);
                obj.modelObj = modelObj;

                % Fetches the default values
                obj.updateMParsArg;
                obj.modelFun = @(mPars, dx)obj.modelObj.reference(mPars,dx);
                obj.modelType = modelObj.modelType;
                obj.dimension = modelObj.dimension;
    %             obj.dimension = length(size(img));
                if strcmp(modelObj.modelType,'discrete')
                    obj.pixelSize = 1;
                end
            end
            obj.loadListener;
        end
        function loadListener(obj)
            addlistener(obj, 'mParsArgModified', @mParsArgModified_callback);
        end
        function updateMParsArg(obj)
            % update the mPars' arguments based on the change of the
            % geometric model
            modelObj = obj.modelObj;
            for k = 1:length(modelObj.parsArgName)
                obj.mPars.(modelObj.parsArgName{k}) = modelObj.(modelObj.parsArgName{k});
            end
        end
        function [ax, img] = plot(obj, mPars, varargin)
            if isempty(mPars)
                if ~isempty(obj.mPars)
                    fn = obj.mPars.name;
                    for k = 1:length(fn)
                        mPars.(fn{k}) = obj.mPars.value(k);
                    end
                end
            end
            p = inputParser;
            p.addParameter('Projection','xyz');
            p.addParameter('roiSize',220);
            p.addParameter('pixelSize',25);
            p.addParameter('sigma',50);
            parse(p, varargin{:});
            results = p.Results;
            projection = results.Projection;
            if isequal(obj.modelType, 'discrete')
                % for discrete models
                fig = figure('Name','discreteModel');
                ax = axes(fig);
                newRefLocs = obj.modelFun(mPars);
                indUpperRef = newRefLocs.z >= 0;
                axp1 = plot3(ax, newRefLocs.x(indUpperRef), newRefLocs.y(indUpperRef), newRefLocs.z(indUpperRef), ' ow', 'MarkerSize', 10, 'MarkerFaceColor', [0.9 0.1 0.1]);
                hold on
                plot3(ax, newRefLocs.x(~indUpperRef), newRefLocs.y(~indUpperRef), newRefLocs.z(~indUpperRef), ' ow', 'MarkerSize', 10, 'MarkerFaceColor', [0.9 0.6 0.6])
                axis(ax,'equal')
                legend(axp1,{'Reference', 'Target'},'Location','northwest')
            else
                % for continuous models
                fig = figure('Name','continuousModel');
                ax = axes(fig);
                pixelSize = results.pixelSize;
                sigma = results.sigma/pixelSize;
                roiSize = results.roiSize;
                bound = roiSize/2;
                rangex = [-bound bound]; % this should be linked to the ROI size
                rangey = rangex;
                rangez = rangex;
                imgDim = ceil((range(rangez))/pixelSize);

                img = voxelblurr(@obj.fun,mPars,sigma, pixelSize,rangex, rangey, rangez);
                img = img./sum(img,1:3);
                if ndims(img)==2
                else
                    switch projection
                        case 'xyz'
                            disp('Projection for 3D model is required.');
                            return
                        case 'xy'
                            img = squeeze(max(img,[],3))';

                        case 'xz'
                            img = squeeze(mean(img,1))';

                        case 'yz'
                            img = squeeze(mean(img,2))';
                    end
                end
                imagesc(ax,img);
                lut = mymakelut('red hot');
                colormap(ax, lut);
%                 caxis(ax, [0 3.5e-6]);
            end
            % move to script
%             indUpper = locs.znm>=0;
%             axp2 = plot3(axs, locs.xnmrot(indUpper), locs.ynmrot(indUpper), locs.znm(indUpper), ' ob', 'MarkerSize', 4, 'MarkerFaceColor', [0.1 0.1 1]);
%             plot3(axs, locs.xnmrot(~indUpper), locs.ynmrot(~indUpper), locs.znm(~indUpper), ' ob', 'MarkerSize', 4, 'MarkerFaceColor', [0.7 0.7 1])
        end
        function [ax, img] = isosurfacePlot(obj, mPars, varargin)
            if isempty(mPars)
                fn = obj.mPars.name;
                for k = 1:length(fn)
                    mPars.(fn{k}) = obj.mPars.value(k);
                end
            end
            p = inputParser;
            p.addParameter('Projection','xyz');
            p.addParameter('roiSize',220);
            p.addParameter('pixelSize',25);
            p.addParameter('sigma',50);
            parse(p, varargin{:});
            results = p.Results;
            
            % get the anchor points of the model
%             sigma = results.sigma/pixelSize;
            img = obj.getImage(mPars,'pixelSize',results.pixelSize);
            patchPlot(img)
            
            img2 = img;
            img2(end/2:end,:,:) = 0;
            img(1:end/2-2,:,:) = 0;
            fig = figure;
            ax = axes(fig);
%             isosurfacePlot(ax, img, 'v2', img2)
            
        end
        
        function [ax, img] = patchPlot(obj, mPars, varargin)
            if isempty(mPars)
                fn = obj.mPars.name;
                for k = 1:length(fn)
                    mPars.(fn{k}) = obj.mPars.value(k);
                end
            end
            p = inputParser;
            p.addParameter('Projection','xyz');
            p.addParameter('roiSize',220);
            p.addParameter('pixelSize',2);
            p.addParameter('sigma',15);
            p.addParameter('lPars',[]);
            p.addParameter('axes',[]);
            p.addParameter('azi_light',45);
            p.addParameter('ele_light',60);
            p.addParameter('azi_view',0);
            p.addParameter('ele_view',-5);
            p.addParameter('colormap','winter');
            p.addParameter('FaceColor','#C0C0C0');
            p.addParameter('isoCutoff',1.1);
            parse(p, varargin{:});
            p = p.Results;
            
            % get the anchor points of the model
%             sigma = results.sigma/pixelSize;
            img = obj.getImage(mPars,'pixelSize',p.pixelSize, 'lPars', p.lPars);
            if isempty(p.axes)
                fig = figure;
                ax = axes(fig);
            else
                ax = p.axes;
            end
            inp = rmfield(p,{'Projection','roiSize','pixelSize','sigma','axes','lPars'});
            inp=[fieldnames(inp).'; struct2cell(inp).'];
            inp = inp(:).';
            patchPlot(ax,img, inp{:})
            
        end
        
        function img = getImage(obj, mPars, varargin)
            saveCache = false;
            mode = 'regular';
            if ~isempty(obj.ParentObject)
                if isempty(obj.modelObj)
                    modelType = '';
                else
                    modelType = obj.modelObj.modelType;
                end
                switch modelType
                    case 'background'
                        mode = 'background';
                        switch obj.ParentObject.status
                            case 'initial'
                                saveCache = true;
                            otherwise
                                saveCache = false;
                        end
                    otherwise
                end
            end
            
            
            % Getting the gaussian filterd image of the geometric model.
            if ~exist('mPars','var')
                mPars = obj.ParentObject.exportPars(obj.ID,'mPar');
            end
            % parse varargin
            p = inputParser;
            p.addParameter('pixelSize', obj.pixelSize);
            p.addParameter('projection', 'none');
            if ~isempty(obj.ParentObject)
                p.addParameter('roiSize', obj.ParentObject.roiSize);
            else
                p.addParameter('roiSize', 300);
            end
            
            if ~isempty(obj.ParentObject)
                p.addParameter('imgExtension', obj.ParentObject.imgExtension);
            else
                p.addParameter('imgExtension', 0);
            end
            p.addParameter('useLocprecnm', ~obj.fixSigma);
            p.addParameter('sigma', []);
            p.addParameter('lPars', []);
            parse(p,varargin{:});
            results = p.Results;
            pixelSize = results.pixelSize;
            sigma = results.sigma;
            roiSize = results.roiSize;
            projection = results.projection;
            imgExtension = results.imgExtension;
            useLocprecnm = results.useLocprecnm;
            
            % Getting the gaussian filterd image
            if useLocprecnm
                sigma = obj.locsPrecFactor/pixelSize;
            else
                if isempty(sigma)
                    sigma = obj.sigma/pixelSize;
                else
                    % then do nothing
                end
            end
            bound = roiSize/2+imgExtension;
            rangex = [-bound bound];
            rangey = rangex;
            rangez = rangex;
            if strcmp(mode,'regular')||saveCache
                imgDim = ceil((range(rangez))/pixelSize);
                if useLocprecnm&&isprop(obj,'pseudoModel')
                    sigma = (obj.pseudoModel.locprecnm+obj.sigmaFactor(2))./pixelSize;
                end
                modelCoord = obj.getPoint(mPars);
                if ~isempty(results.lPars)
                    modelCoord = obj.ParentObject.transModPoint(modelCoord,'lPars',results.lPars);
                end
                switch projection
                    case 'none'
                    case 'xy'
                        modelCoord = rmfield(modelCoord,'z');
                    case 'xz'
                        modelCoord = obj.getPoint(mPars);
                        modelCoord.y = modelCoord.z;
                        modelCoord = rmfield(modelCoord,'z');
                    case 'yz'
                        modelCoord = obj.getPoint(mPars);
                        modelCoord.x = modelCoord.z;
                        modelCoord = rmfield(modelCoord,'z');
                end
                img = voxelblurr_yw(modelCoord,mPars,sigma, pixelSize,rangex, rangey, rangez,obj.modelType);
                img = img./sum(img,1:3);
                if saveCache
                    obj.img = img;
                end
            else
                img = obj.img;
            end
            if isequal(mode,'background')
                if 1
                    oriVal = obj.ParentObject.allParsArg.value;
                    obj.ParentObject.allParsArg.value(~obj.ParentObject.allParsArg.fix) = obj.ParentObject.currentFitPars;
                    mPar_m1= obj.ParentObject.exportPars(1,'mPar');
    %                 lPar_m1= obj.ParentObject.exportPars(1,'lPar');
                    obj.ParentObject.allParsArg.value = oriVal;
                    m1Points = obj.ParentObject.model{1}.getPoint(mPar_m1);
                    m1Points.x = m1Points.x-obj.pixelSize;
                    m1Points.y = m1Points.y-obj.pixelSize;
                    m1Points.z = m1Points.z-obj.pixelSize;
    %                 m1Locs.xnm = m1Points.x;
    %                 m1Locs.ynm = m1Points.y;
    %                 m1Locs.znm = m1Points.z;
    %                 m1Locs = obj.ParentObject.locsHandler(m1Locs,lPar_m1,[],'usedformalism','rotationMatrixRev');
    %                 m1Points.x = m1Locs.xnm;
    %                 m1Points.y = m1Locs.ynm;
    %                 m1Points.z = m1Locs.znm;

    %                 mxyz = [-(obj.ParentObject.roiSize)/2 (obj.ParentObject.roiSize+1)/2];
                    idxSet2Zero = voxelblurr_yw(m1Points,mPars,sigma, pixelSize,rangex, rangey, rangez,obj.modelType);
                    temp = idxSet2Zero(idxSet2Zero>0);
                    cutoff = prctile(temp(:),75);
                    idxSet2Zero = idxSet2Zero>cutoff;
                    img(idxSet2Zero) = 0;
                end
                img = img./sum(img,1:3);
            end
        end
        function [ref,mp] = getPoint(obj, mPars, varargin)
            % Getting sampled points from the model.
            p = inputParser;
            % p.addParameter('factor', single(min(obj.sigma)./obj.pixelSize));
            % factor 0.75 means the sampled points will be 0.75*sigma away
            % from each other.
            p.addParameter('factor', []);
            parse(p,varargin{:});
            results = p.Results;
            
            % dx, the lower the denser, is a parameter controlling the
            % sampling rate. dx means the minimal distance between nearst
            % sampled points.
            if isempty(results.factor)
                if obj.fixSigma
                     dx = obj.ParentObject.refPoint_spacing;
    %                 dx = obj.refPoint_spacing;
                else
                    if strcmp(obj.modelType, 'continuous')
                        %!!! This might be wrong. dx should not be dependent on
                        % obj.locsPrecFactor.
                        dx = obj.ParentObject.refPoint_spacing;
    %                     dx = obj.ParentObject.refPoint_spacing;
                    else
                        if isempty(obj.ParentObject)
                            dx = obj.samplingFactor;
                        else
                            dx = obj.ParentObject.refPoint_spacing;
                        end
                    end
                end
            else
                dx = results.factor;
            end
            if obj.dimension == 3
                [ref.x,ref.y,ref.z,ref.n, mp] = obj.fun(mPars,dx);
                ref.z = ref.z(ref.n~=0);
            else
                [ref.x,ref.y,~,ref.n, mp] = obj.fun(mPars,dx);
            end
                ref.x = ref.x(ref.n~=0);
                ref.y = ref.y(ref.n~=0);
                ref.n = ref.n(ref.n~=0);
        end
        function [x,y,z,n,mp] = fun(obj, mPars, dx) % convert the output of model for voxelblurr          
            [model,mp] = obj.modelFun(mPars, dx);
            x = model.x;
            y = model.y;
            if isfield(model, 'z')
                z = model.z;
            else
                z = [];
            end
            n = model.n;
        end

        function rmMPars(obj,mPar)
            obj.mPar.name
            obj.mPars.name = modelObj.name;
            obj.mPars.fix = modelObj.fix;
            obj.mPars.value = modelObj.value;
            obj.mPars.lb = modelObj.lb;
            obj.mPars.ub = modelObj.ub;
            obj.mPars.min = modelObj.min;
            obj.mPars.max = modelObj.max;
        end
        function mParsArgModified_callback(obj,b)
            obj.updateMParsArg
        end
        
        function [sigmaFactor, sigmaSet, sigmaZSet] = deriveSigma(obj, locs, varargin)
            % :meth:`deriveSigma` derives the final sigma used for
            % fitting. When :attr:`fixSigma` is set as true, sigma are
            % derived based on pre-defined values. Otherwise, sigma are
            % derived based on localization precisions. For a continuous
            % model, the minimum sigma is defined as the median of
            % localization precisions.
            % 
            % Uasage:
            %   obj.deriveSigma(locs)
            %
            % Args:
            %   obj (functionModel object): an object created by
            %   :func:`functionModel`.
            %   locs (structure array): a typical localization structure
            %   array used in SMAP.
            %
            % Returns:
            %   sigmaFactor (numeric vector): a 1-by-2 vector that
            %   determines the fold of localization precisions used for
            %   fitting.
            %   sigmaSet (numeric vector | numeric scalar): sigma used for
            %   fitting. A N-by-1 vector, where N is the number of
            %   localiztions when :attr:`fixSigma` is true.
            %   sigmaZSet (numeric vector | numeric scalar): z sigma used for
            %   fitting. A N-by-1 vector, where N is the number of
            %   localiztions when :attr:`fixSigma` is true.
            %
            % Last update:
            %   28.04.2022
            %
            % See also:
            %   :class:`functionModel`

            inp = inputParser;
            minSigma = obj.ParentObject.getAdvanceSetting('minSigma');
            if isempty(minSigma)
                inp.addParameter('minSigma', 'median')
            else
                inp.addParameter('minSigma', minSigma)
            end
            inp.parse(varargin{:})
            inp = inp.Results;

            if ~obj.fixSigma
                sigmaFactor = obj.sigmaFactor;
                
                %% sigmaFactor contains 2 elements. The first one is the factor itself, the second is an offset.
                % By default the offset is 0 unless specified otherwise.
                if length(sigmaFactor)==1
                    sigmaFactor(1,2) = 0;
                end
                sigmaSet = locs.locprecnm;
                if obj.dimension == 3
                    sigmaZSet = locs.locprecznm;
                end

                switch obj.modelType
                    case 'discrete'
                        % do nothing here if discrete
                    otherwise
                        if strcmp(inp.minSigma,'median')
                            % define the median locprec as the minimum
                            medSig = median(locs.locprecnm);
                            sigmaSet(sigmaSet<medSig)=medSig;
                            if obj.dimension == 3
                                medSigz = median(locs.locprecznm);
                                sigmaZSet(sigmaZSet<medSigz) = medSigz;
                            end
                        elseif strcmp(inp.minSigma,'off')
                            % do nothing here if off
                        end
                end
            else
                sigmaFactor = [1 0];
                sigmaSet = obj.sigma;
                if obj.dimension == 3
                    sigmaZSet = sigmaSet*obj.sigmaZFactor;
                end
            end
            obj.sigmaFactor = sigmaFactor;
            obj.sigmaSet = sigmaSet;
            if obj.dimension == 3
                obj.sigmaZSet = sigmaZSet;
            end
        end

        %% define the dependent parameters
        function locsPrecFactor = get.locsPrecFactor(obj)
            % Define the locsPrecFactor based on either obj.sigmaSet or
            % obj.sigma. obj.sigmaSet has a higher priority.
            if isempty(obj.sigmaSet)
                locsPrecFactor = obj.sigma;
            else
                locsPrecFactor = min(obj.sigmaSet);
            end
        end
        
    end
    methods(Access = protected)
        function cp = copyElement(obj)
            cp = copyElement@matlab.mixin.Copyable(obj);
            cp.modelObj = copy(cp.modelObj);
            cp.modelObj.ParentObject = cp;
        end
    end
end
