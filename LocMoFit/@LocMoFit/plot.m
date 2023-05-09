function [ax,finalImg] = plot(obj,locs,varargin)
% :meth:`plot` visualize the fit result.
%
% Args:
%   obj (:obj:`LocMoFit`): an object created by :meth:`LocMoFit`.
%   locs (structure array): the typical localization structure array used
%   in SMAP.
%
% Name-Value Arguments:
%   bestPar
%   lPars
%   mPars
%   plotType
%   whichModel
%   Projection
%   pixelSize (numeric scalar): the pixel size of the rendered image. Default: determined by the :obj:`SMLMModel`.
%   axes
%   movModel
%   normalizationMax
%   shift
%   color
%   doNotPlot
%   modelSamplingFactor
%   displayLocs
%   color_max
%   MarkerFaceAlpha_Locs (numeric scalar): this determines the transparency of the locs points. Default: 0.5.
%
% Returns:
%   ax (Axes object): the axes object where the fit is visualized.
%   finalImg (structure array, numeric array): either a rendered image (numeric array) or a structure array representing the current model.
%   
% Last update:
%   06.05.2022
%

%% PLOT Display the fit results
% For a coordinate-based model in 3D, the fitted model can be displayed as points, 
% projected images, or both simultaneously.

%% Initiation
p = inputParser;
p.addParameter('bestPar',true);
p.addParameter('lPars',[]);
p.addParameter('mPars',[]);
p.addParameter('plotType','image');                   % can be either image or point
p.addParameter('whichModel', obj.numOfModel);
p.addParameter('Projection', 'xy');
p.addParameter('pixelSize', obj.model{1}.pixelSize);
p.addParameter('axes', []);
p.addParameter('movModel', false);
p.addParameter('normalizationMax', true);
p.addParameter('shift', [0 0 0]);                     %!!! This is implemented for the case when model is cropped at the edge. Need to be further tested.
p.addParameter('color', {});                          %!!! need to be further introduced. take the input from SMAP
p.addParameter('doNotPlot', false);
p.addParameter('modelSamplingFactor', []);
p.addParameter('displayLocs',1)
p.addParameter('alpha_Locs', 0.2);
p.addParameter('color_max',ones([1 obj.numOfLayer])*255)
% If the first argument is an handle of axes, then use it as the parent
if ~exist('locs',"var")
    locs = obj.locs;
end
if isa(locs,'matlab.graphics.axis.Axes')
    ax = locs;
    locs = varargin{1};
    varargin(1)=[];
end
parse(p, varargin{:});
results = p.Results;
whichModel = results.whichModel;
projection = results.Projection;
pixelSize = results.pixelSize;
%% Check the model types and the number of layers
% If any of the model is an image then display an image at the end.
for k = obj.numOfModel:-1:1
    modelType{k} = obj.model{k}.modelType;
    modelLayer(k) = obj.model{k}.layer;
end
if ismember({'image'}, modelType)
    results.plotType = 'image';
end
allModelLayers = unique(modelLayer);
%% Transform lPars to the shared coordinates system
% Only do lPars have additionality always move the model rather than the locs
if results.bestPar
    lPars = {};
    for k = 1:obj.numOfModel
        indModel = obj.allParsArg.model == k;
        indLp = ismember(obj.allParsArg.type,'lPar');
        indMp = ismember(obj.allParsArg.type,'mPar');
        % lPars
        fn = obj.allParsArg.name(indModel&indLp);
        bestFit = obj.allParsArg.value(indModel&indLp);
%         if k == 1
            for l = 1:length(fn)
                lPars{k}.(fn{l}) = bestFit(l);
            end
%         else
            % aligned all models to the model one
%             for l = 1:length(fn)
%                 lPars{k}.(fn{l}) = lPars{1}.(fn{l}) + bestFit(l);
%             endR
%         end
    end
else
    lPars = results.lPars;
end
%% Get images of the models
% Image model should always display image visulization, while point model can 
% display either image or point visulization.
% Determine the image size
modelType = obj.model{1}.modelType;
modelDim = obj.model{1}.dimension;
imgSize = size(obj.model{1}.img);
if ~strcmp(modelType,'image')
    imgSize = imgSize./pixelSize;
end
lParM1 = obj.exportPars(1,'lPar');
if strcmp(modelType,'image')
    lParM1.x = lParM1.x-pixelSize;
    lParM1.y = lParM1.y-pixelSize;
    if obj.dataDim==3
        lParM1.z = lParM1.z-pixelSize;
    end
end
if ~isempty(locs)
    % do not transform the locs if moving the model
    if ~results.movModel
        newLocs = obj.locsHandler(locs, lParM1,1);
    else
        newLocs = locs;
    end
end
if any(imgSize == 0)
    imgSize = repelem((obj.roiSize+2*obj.imgExtension), modelDim)./pixelSize;
end

%% Get each model
%         For the image type
if isequal(results.plotType,'image')
    % for an image model, use the method 'plot()' of image model class
    % to generate the image
    artBoard = zeros(imgSize);
    for k = 1:obj.numOfModel
        if ~isequal(obj.model{k}.modelType,'image')
            obj.model{k}.deriveSigma(locs);
        end
        if k==1
            oneMPars = obj.exportPars(k,'mPar');
            if ~results.movModel
                modelImage{k} = obj.model{k}.getImage(oneMPars,'pixelSize', pixelSize,'roiSize',obj.roiSize);
            else %!!! only for the tac figures
                oneLPars = obj.exportPars(k,'lPar');
                img = obj.model{k}.getImage(oneMPars,'pixelSize', pixelSize,'roiSize',obj.roiSize);
                [X,Y] = meshgrid(1:imgSize(1),1:imgSize(2));
                F = griddedInterpolant(img,'cubic', 'nearest');
                xPos=oneLPars.x;
                yPos=oneLPars.y;
                zrot=oneLPars.zrot;
                XMean = mean(X(:)); YMean = mean(Y(:));
                cenX = X-XMean;
                cenY = Y-YMean;
                cenX = cenX-(xPos)./pixelSize;
                cenY = cenY-(yPos)./pixelSize;
                [cenX, cenY] = rotcoord(cenX(:),cenY(:),-zrot*pi/180);
                valNew = F(cenX+XMean + results.shift(2)/pixelSize , cenY+YMean+results.shift(1)/pixelSize);
                ind = sub2ind(imgSize, X(:),Y(:));
                finalImg = artBoard;
                finalImg(ind)=valNew;
                modelImage{k} = finalImg;
                
            end
        else
            oneMPars = obj.exportPars(k,'mPar');
            img = obj.model{k}.getImage(oneMPars,'pixelSize', pixelSize,'roiSize',obj.roiSize);
            
            % translate lPar to mPar
            oneLPars = obj.exportPars(k,'lPar');          % get lPars
            
            % get unit grid points
            if isfield(newLocs,'znm')
                % for 3D
                [X,Y,Z] = meshgrid(1:imgSize(1),1:imgSize(2),1:imgSize(3));
                zPos=oneLPars.z;
                xrot=oneLPars.xrot;
                yrot=oneLPars.yrot;
                
                ZMean = mean(Z(:));
                cenZ = Z-ZMean;
            else
                % for 2D
                [X,Y] = meshgrid(1:imgSize(1),1:imgSize(2));
            end
            F = griddedInterpolant(img,'cubic', 'nearest');
            
            xPos=oneLPars.x;
            yPos=oneLPars.y;
            zrot=oneLPars.zrot;
            
            XMean = mean(X(:)); YMean = mean(Y(:));
            cenX = X-XMean;
            cenY = Y-YMean;
            
            
            if isfield(newLocs,'znm')
                % for 3D
                [cenX, cenY, cenZ] = rotcoord3(cenX(:),cenY(:),cenZ(:),deg2rad(xrot),deg2rad(yrot),deg2rad(zrot),'XYZ');
                valNew = F(cenX+XMean-(xPos)./pixelSize, cenY+YMean-(yPos)./pixelSize, cenZ+ZMean-(zPos)./pixelSize);
                ind = sub2ind(imgSize, X(:),Y(:),Z(:));
            else
                % for 2D
                [cenX, cenY] = rotcoord(cenX(:),cenY(:),-zrot*pi/180);
                valNew = F(cenX+XMean-(xPos)./pixelSize + results.shift(2)/pixelSize , cenY+YMean-(yPos)./pixelSize+results.shift(1)/pixelSize);
                ind = sub2ind(imgSize, X(:),Y(:));
            end
            finalImg = artBoard;
            finalImg(ind) = valNew;
            modelImage{k} = finalImg;
        end
    end
else
    % For the point type
%     modPoint = obj.getModPoint(results.modelSamplingFactor);
end
%% Render the models
% From models to layers: generate for each layer an image
dataCol = {'r','g','b'};
if isequal(results.plotType,'image')
%     nameAllLut = mymakelut;
    img2disp = artBoard(:,:,1);
    for ch = 1:length(allModelLayers)
        oneCh = artBoard;
        indModelOneCh = find(modelLayer == allModelLayers(ch));
        if isempty(obj.linkedGUI)
            warning('The GUI object is not linked.')
            nameLut{ch} = 'cyan';
        else
            if obj.numOfLayer==1
                nameLut{ch} = 'cyan';
            else
                nameLut{ch} = obj.linkedGUI.getPar(['layer' num2str(ch) '_lut']).selection;
            end
        end
        for k = 1:length(indModelOneCh)
            indOneModel = indModelOneCh(k);
            if ~isempty(obj.fitInfo)
                % before 200605
%                 oneImage = modelImage{indOneModel}.*obj.fitInfo.weightModel{indOneModel};
%                 disp('!!!200605: changed, check here if warning.')
                oneImage = modelImage{indOneModel};
            else
                oneImage = modelImage{indOneModel};
            end
            if k == 1
                oneCh = artBoard;
            end
            oneCh = oneCh+oneImage;
        end
        if obj.dataDim == 3
            switch projection
                case 'xy'
                    oneCh = squeeze(mean(oneCh,3));
                case 'xz'
                    oneCh = squeeze(mean(oneCh,2));
                case 'yz'
                    oneCh = squeeze(mean(oneCh,1));
            end
        end
        if ~isempty(obj.weightLayer)
            % before 200605
%               img2disp = img2disp + ind2rgb(ceil((oneCh./max(oneCh,[],1:length(size(oneCh)))).*obj.weightLayer(ch)*255)', mymakelut(nameAllLut{ch})); %%% !!!!!!quick and dirty
%               disp('!!!200605: changed, check here if warning.')
            img2disp = img2disp + ind2rgb(ceil((oneCh./max(oneCh,[],1:length(size(oneCh)))).*results.color_max(ch))', mymakelut(nameLut{ch})); %%% !!!!!!quick and dirty
        else
            img2disp = img2disp + ind2rgb(ceil((oneCh./max(oneCh,[],1:length(size(oneCh)))).*results.color_max(ch))', mymakelut(nameLut{ch}));
        end
    end
    finalImg = img2disp; %%%quick and dirty
else
    if isempty(obj.modelLayer)
        obj.updateLayer
    end
    layerPoint = obj.getLayerPoint(results.modelSamplingFactor);
end
%% Display images
if isequal(results.plotType,'image')&&~isempty(locs)
    % project the 3D image
    if obj.dataDim == 3
        switch projection
            case 'xy'
                dim1 = newLocs.xnm./pixelSize+imgSize(2)/2;
                dim2 = newLocs.ynm./pixelSize+imgSize(1)/2;
            case 'xz'
                dim1 = newLocs.xnm./pixelSize+imgSize(2)/2;
                dim2 = newLocs.znm./pixelSize+imgSize(3)/2;
            case 'yz'
                dim1 = newLocs.ynm./pixelSize+imgSize(1)/2;
                dim2 = newLocs.znm./pixelSize+imgSize(3)/2;
        end
    end
    % show the image
    if ~exist('ax', 'var')
        if ~results.doNotPlot
            fig = figure;
            ax = axes(fig);
        else
            ax = [];
        end
    end
    if ~results.doNotPlot
        imagesc(ax, img2disp);
    end
    
    % display the locs channel by channel
    if results.displayLocs
        if ~results.doNotPlot
            for layer = 1:length(allModelLayers)
                % get the locs too
                hold(ax, 'on')
                if ~isempty(locs) && obj.dataDim==2
                    dim1 = newLocs.xnm./pixelSize+imgSize(2)/2;
                    dim2 = newLocs.ynm./pixelSize+imgSize(1)/2;
                end
                currentLut = mymakelut(nameLut{layer});
                if obj.numOfLayer==1
                    oneLut = mymakelut('red hot');
                    oneColor = oneLut(140,:);
                    plot(ax, dim1(locs.layer==layer), dim2(locs.layer==layer), ' o', 'MarkerFaceColor', oneColor, 'MarkerEdgeColor','k', 'MarkerSize',2.5);
                else
                    plot(ax, dim1(locs.layer==layer), dim2(locs.layer==layer), ' o', 'MarkerFaceColor', currentLut(100,:), 'MarkerEdgeColor','w', 'LineWidth',1, 'MarkerSize',3);
                end
                hold(ax,'off')
            end
        end
    end
else
    % for point view
    if ~exist('ax', 'var')
        if ~results.doNotPlot
            fig = figure;
            ax = axes(fig);
        else
            ax = [];
        end
    end
    if ~isfield(locs,'znm')&&~isempty(locs)
        newLocs.znm = zeros(size(newLocs.xnm));
    end
    for layer = 1:length(allModelLayers)
        if ~isfield(layerPoint{layer},'z')
            
        elseif isempty(layerPoint{layer}.z )
            disp('Check here when getting error.')
            % 200610 modified:
%             layerPoint{layer}.z = zeros(size(layerPoint{layer}.x));
            layerPoint{layer} = rmfield(layerPoint{layer}, 'z');
        end
        if ~results.doNotPlot
            hold(ax, 'on')
            scatter3(ax,layerPoint{layer}.x,layerPoint{layer}.y,layerPoint{layer}.z, 16, 'filled', 'o', 'MarkerFaceColor', dataCol{layer}, 'MarkerEdgeColor','k', 'LineWidth',1)
            scatter3(ax,newLocs.xnm(locs.layer==layer),newLocs.ynm(locs.layer==layer),newLocs.znm(locs.layer==layer),16, 'filled', 'o', 'MarkerFaceColor', dataCol{layer}, 'MarkerEdgeColor','w', 'LineWidth',0.5, 'MarkerFaceAlpha', results.alpha_Locs, 'MarkerEdgeAlpha', results.alpha_Locs)
            hold(ax, 'off')
        end
    end
    % display all others
    %         newLocs.xnm(~ismember(locs.layer, allModelLayers))
    finalImg = layerPoint;
end
if ~results.doNotPlot
    axis(ax,'equal')
end
end