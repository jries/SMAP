function intensityVal = modelHandler(obj, locs, mPars, varargin)
%% MODELHANDLER Get intensity of the model at locs' coordinates given a set of model parameters
p = inputParser;
p.addParameter('lPars',[]);
if ~isempty(obj.ParentObject)
    p.addParameter('gaussDistMode',obj.ParentObject.getAdvanceSetting('gaussDistMode'));
    p.addParameter('gaussDistCutoff',obj.ParentObject.getAdvanceSetting('gaussDistCutoff'));
else
    p.addParameter('gaussDistMode','ordinary');
    p.addParameter('gaussDistCutoff',3.5);
end
parse(p,varargin{:});
lPars = p.Results.lPars;
p = p.Results;

%% here for setting the sigma
switch obj.dimension
    case 2
        [sigmaFactor, sigmaSet] = obj.deriveSigma(locs);
    case 3
        [sigmaFactor, sigmaSet, sigmaZSet] = obj.deriveSigma(locs);
end
    
if ~any([isequal(obj.modelType, 'continuous') isequal(obj.modelType, 'background')])
    %% For a discrete model
    %         refLocs = obj.modelFun(mPars);
    
    % This term determines the sampling rate of the model
    refLocs = obj.getPoint(mPars);
    
    lParRefSeries = unique([1 obj.ID]);
    if obj.dimension == 3
        % sequentially apply the transformation
        refLocs.xnm = refLocs.x; refLocs.ynm = refLocs.y; refLocs.znm = refLocs.z;
       
        for k = 1:length(lParRefSeries)
            if obj.ID == 1 && ~isempty(obj.ParentObject.getTemp('freeRot')) && obj.ParentObject.getTemp('freeRot')
                usedformalism = 'lieAlgebra';
            else
                usedformalism = 'rotationMatrixRev';
            end
            refLocs = obj.ParentObject.locsHandler(refLocs, lPars{lParRefSeries(k)}, 0, 'order_transform', 'RT', 'usedformalism', usedformalism);
        end
        refLocs.x = refLocs.xnm; refLocs.y = refLocs.ynm; refLocs.z = refLocs.znm;
    else
        for k = 1:length(lParRefSeries)
            zrot = deg2rad(lPars{lParRefSeries(k)}.zrot);
            [x,y] = rotcoord2(refLocs.x,refLocs.y,zrot');
            x = x';
            y = y';
            refLocs.x = x'+lPars{lParRefSeries(k)}.x;
            refLocs.y = y'+lPars{lParRefSeries(k)}.y;
            refLocs.n = repelem(refLocs.n,1, size(lPars{lParRefSeries(k)}.y,2));
        end
%         lInRoi = refLocs.x.^2+refLocs.y.^2 <=(obj.ParentObject.roiSize/2)^2;
%         refLocs.x = refLocs.x(lInRoi); refLocs.y = refLocs.y(lInRoi); refLocs.n = refLocs.n(lInRoi);
    end
    
    % Deal with the gaussDistCutoff for the fast-mode gaussDist
    fitter = obj.ParentObject;
    oldDistCutoff = fitter.getTemp('oldDistCutoff');
    if isempty(oldDistCutoff)
        fitter.setTemp('oldDistCutoff', p.gaussDistCutoff)
        fitter.setTemp('sumGaussIntensity', 1-(1-mvncdf(repelem(((p.gaussDistCutoff^2)/3)^(1/2),1,3),[0 0 0],diag([1 1 1])))*2)
    else
        if oldDistCutoff~=p.gaussDistCutoff
            fitter.setTemp('oldDistCutoff', p.gaussDistCutoff)
            fitter.setTemp('sumGaussIntensity', 1-(1-mvncdf(repelem(((p.gaussDistCutoff^2)/3)^(1/2),1,3),[0 0 0],diag([1 1 1])))*2)
        end
    end
    sumGaussIntensity = fitter.getTemp('sumGaussIntensity');
    
    if obj.dimension == 3
        % for 3D
        intensityVal = gaussDist(refLocs,locs.ynm, locs.xnm, locs.znm, sigmaSet+sigmaFactor(2),sigmaSet+sigmaFactor(2),sigmaZSet, 'sigFactor', sigmaFactor(1), 'distMode',p.gaussDistMode,'distCutoff',p.gaussDistCutoff)./sum(refLocs.n,1);
        likelihoodBoun = 0;
    else
        % for 2D
        intensityVal = gaussDist2(refLocs,locs.ynm, locs.xnm, sigmaSet+sigmaFactor(2), sigmaSet+sigmaFactor(2), 'sigFactor', sigmaFactor(1))./sum(refLocs.n,1);
        likelihoodBoun = 0;
    end
    %                 likelihoodBoun = gaussDist(refLocs,obj.sigma*2,obj.sigma*2,obj.sigma*2, obj.sigma,obj.sigma,obj.sigma);
    intensityVal(intensityVal<likelihoodBoun) = 0;
else
    roiSize = obj.ParentObject.roiSize;
    imgExtension = obj.ParentObject.imgExtension;
    sigma = obj.locsPrecFactor/obj.pixelSize; % was 12/pixelSize
    bound = roiSize/2+imgExtension;
    rangex = [-bound bound];
    rangey = rangex;
    rangez = rangex;
    %% For a continous model
    % Generate images given parameters
    % mimic vectorized function
    imgDim = ceil((range(rangez))/obj.pixelSize);
    if isempty(mPars)
        numOfParSet = 1;
    else
        fn = fieldnames(mPars);

        for l = length(fn):-1:1
            numOfParSet(l) = length(mPars.(fn{l}));
        end
    end
    
    if obj.dimension == 3
        img = zeros([imgDim imgDim imgDim max(numOfParSet)]);
    else
        img = zeros([imgDim imgDim max(numOfParSet)]);
    end
    % create a stack of image
    if any(numOfParSet>1)
        for k = 1:length(mPars.(fn{1}))     % the kth set of parameters
            for l = 1:length(fn)
                if length(mPars.(fn{l}))>1
                    oneMPars.(fn{l}) = mPars.(fn{l})(1,k);
                else
                    oneMPars.(fn{l}) = mPars.(fn{l});
                end
            end
            oneImg = obj.getImage(oneMPars);
            if ndims(oneImg) == 3
                img(:,:,:,k) = oneImg;
            else
                img(:,:,k) = oneImg;
            end
            
        end
    else
        oneImg = obj.getImage();
        if ndims(oneImg) == 3
            img(:,:,:) = oneImg;
        else
            img(:,:) = oneImg;
        end
    end
    % normalization of the image
    img = (img./sum(img,1:obj.dimension))./obj.pixelSize^obj.dimension;                  % normalize the images
    if ndims(oneImg) == 3
        %             img = img./sum(img,1:3).*(1-obj.ParentObject.roiSize^2*obj.eps);
        zeroInd = abs(locs.xnm)>bound|abs(locs.ynm)>bound|abs(locs.znm)>bound;
        layerInd = locs.layer ~= obj.layer;                                     % find the locs that are not in the layer of the model
        zeroInd=layerInd|zeroInd;
    else
        %             img = img./sum(img,1:2).*(1-obj.ParentObject.roiSize^2*obj.eps);
        zeroInd = abs(locs.xnm)>bound|abs(locs.ynm)>bound;
        layerInd = locs.layer ~= obj.layer;
        zeroInd=layerInd|zeroInd;
    end
    % Get intensity at locs' coordinates
    intensityVal = zeros(size(locs.xnm));
    for k = 1:size(locs.xnm,2)
        if size(img, obj.dimension+1)==1
            l = 1;
            %% !!!not yet finished
            if ndims(oneImg) == 3
                % for 3D
                if sum(~zeroInd(:,k))>0
                    intensityVal(~zeroInd(:,k),k) = ba_interp3(double(squeeze(img(:,:,:,l))),...    % x and y need to be exchange to match the indexing of a matrix
                        double((locs.ynm(~zeroInd(:,k),k)+bound)/obj.pixelSize),...
                        double((locs.xnm(~zeroInd(:,k),k)+bound)/obj.pixelSize),...
                        double((locs.znm(~zeroInd(:,k),k)+bound)/obj.pixelSize),'cubic');
                end
            else
                % for 2D
                if sum(~zeroInd(:,k))>0
                    intensityVal(~zeroInd(:,k),k) = interp2(double(squeeze(img(:,:,l))),...         % x and y need to be exchange to match the indexing of a matrix
                        double((locs.ynm(~zeroInd(:,k),k)+bound)/obj.pixelSize),...
                        double((locs.xnm(~zeroInd(:,k),k)+bound)/obj.pixelSize),'cubic');
                end
            end
        else
            if ndims(oneImg) == 3
                % for 3D
                if sum(~zeroInd(:,k))>0
                    intensityVal(~zeroInd(:,k),k) = ba_interp3(double(squeeze(img(:,:,:,k))),...    % x and y need to be exchange to match the indexing of a matrix
                        double((locs.ynm(~zeroInd(:,k),k)+bound)/obj.pixelSize),...
                        double((locs.xnm(~zeroInd(:,k),k)+bound)/obj.pixelSize),...
                        double((locs.znm(~zeroInd(:,k),k)+bound)/obj.pixelSize),'cubic');
                end
            else
                % for 2D
                if sum(~zeroInd(:,k))>0
                    intensityVal(~zeroInd(:,k),k) = interp2(double(squeeze(img(:,:,k))),...         % x and y need to be exchange to match the indexing of a matrix
                        double((locs.ynm(~zeroInd(:,k),k)+bound)/obj.pixelSize),...
                        double((locs.xnm(~zeroInd(:,k),k)+bound)/obj.pixelSize),'cubic');
                end
            end
        end
    end
    if 0
        figure; imagesc(squeeze(mean(img,3)))
        hold on
        plot(gca, double((locs.ynm(~zeroInd(:,k),k)+bound)/obj.pixelSize), double((locs.xnm(~zeroInd(:,k),k)+bound)/obj.pixelSize),' o')
    end
    intensityVal(zeroInd) = 0;
    %         refLocs.x = 0; refLocs.y = 0; refLocs.z = 0;
    %                 likelihoodBoun = (1/(roiSize/pixelSize)^3)*gaussDist(refLocs,sigma*2,sigma*2,sigma*3*2, sigma,sigma,sigma*3);
    %                 intensityVal(intensityVal<likelihoodBoun)=0;
    intensityVal(isnan(intensityVal)|intensityVal<0)=0;
    %                 intensityVal = intensityVal.*(sum(intensityVal>0)./size(intensityVal,1))+1e-20;
end

