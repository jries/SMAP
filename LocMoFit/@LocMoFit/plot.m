function [ax,finalImg] = plot(obj,locs,varargin)          % visualize the best fit
    % initiation
    p = inputParser;
    p.addParameter('bestPar',true);
    p.addParameter('lPars',[]);
    p.addParameter('mPars',[]);
    p.addParameter('plotType','image');
    p.addParameter('whichModel', obj.numOfModel);
    p.addParameter('Projection', 'xy');
    p.addParameter('pixelSize', 2);
    p.addParameter('axes', []);
    p.addParameter('movModel', true);
    p.addParameter('normalizationMax', true);
    parse(p, varargin{:});
    results = p.Results;
    whichModel = results.whichModel;
    projection = results.Projection;
    pixelSize = results.pixelSize;
    
    
    if results.bestPar                                  %??? maybe user would like to view with specific parameters?
         %% transformed lPars to the shared coordinates system
         % only do lPars have additionality 
         % always move the model rather than the locs
         lPars = {};
         for k = 1:obj.numOfModel
            indModel = obj.allParsArg.model == k;
            indLp = ismember(obj.allParsArg.type,'lPar');
            indMp = ismember(obj.allParsArg.type,'mPar');
            % lPars
            fn = obj.allParsArg.name(indModel&indLp);
            bestFit = obj.allParsArg.value(indModel&indLp);
            if k == 1
                for l = 1:length(fn)
                    lPars{k}.(fn{l}) = bestFit(l);
                end
            else
                for l = 1:length(fn)
                    lPars{k}.(fn{l}) = lPars{k-1}.(fn{l}) + bestFit(l);
                end
            end
         end
    else
        lPars = results.lPars;
    end
        %% get an image of the model
        modelType = obj.model{1}.modelType;
        modelDim = obj.model{1}.dimension;
        imgSize = size(obj.model{1}.img);
        imgSize = imgSize./pixelSize;
        if any(imgSize == 0)
            imgSize = repelem(obj.roiSize, modelDim)./pixelSize;
        end
        finalImg = zeros(imgSize);

        if isequal(modelType, 'discrete')
            fig = figure('Name', 'Fit result');
            ax = axes(fig);
            % need to allow different model
            mPars = obj.exportPars(1,'mPar');
            oneModel = obj.model{1}.modelFun(mPars);
            % translate lPars to mPars
            lPars = obj.exportPars(1,'lPar');
            if results.movModel
                locs = obj.locsHandler(locs, lPars);
                %% !!!specific for NPC, should be removed later
                lLocsUpRing = locs.znm>=0;
                lModelUpRing = oneModel.z >= 0;
                plot3(ax, oneModel.x(lModelUpRing), oneModel.y(lModelUpRing), oneModel.z(lModelUpRing), ' or', 'MarkerSize', 12, 'MarkerFaceColor', [1 0 0]);
                hold(ax, 'on')
                plot3(ax, oneModel.x(~lModelUpRing), oneModel.y(~lModelUpRing), oneModel.z(~lModelUpRing), ' or', 'MarkerSize', 12, 'MarkerFaceColor', [1 0.7 0.7]);
                hold(ax, 'off')
            else
                oneModel.x = oneModel.x + lPars.x;
                oneModel.y = oneModel.y + lPars.y;
                oneModel.z = oneModel.z + lPars.z;
                [oneModel.x, oneModel.z] = rotcoord(oneModel.x, oneModel.z, lPars.yrot*pi/180);
                [oneModel.x, oneModel.y] = rotcoord(oneModel.x, oneModel.y, lPars.zrot*pi/180);
                [oneModel.y, oneModel.z] = rotcoord(oneModel.y, oneModel.z, lPars.xrot*pi/180);
                plot3(ax, oneModel.x, oneModel.y, oneModel.z, ' or')
            end
        elseif isequal(results.plotType, 'point')
            if modelDim == 3
                for k = 1:obj.numOfModel
                    if results.bestPar
                        oneMPars = obj.exportPars(k,'mPar');          % get mPars
                    else
                        oneMPars = results.mPars{k};
                    end
                    reference{k} = obj.model{k}.modelFun(oneMPars,3);

                    if results.movModel
                        locs = obj.locsHandler(locs, lPars{k});
                    else
                        % translate lPars to mPars
                        xPos=lPars{k}.x;
                        yPos=lPars{k}.y;
                        zPos=lPars{k}.z;
                        xrot=lPars{k}.xrot;
                        yrot=lPars{k}.yrot;
                        zrot=lPars{k}.zrot;

                        [reference{k}.x,reference{k}.z] = rotcoord(reference{k}.x+xPos,reference{k}.z+zPos,yrot*pi/180);
                        [reference{k}.y,reference{k}.z] = rotcoord(reference{k}.y+yPos,reference{k}.z,xrot*pi/180);
                        [reference{k}.x,reference{k}.y] = rotcoord(reference{k}.x,reference{k}.y,zrot*pi/180);
                        reference{k}.x = reference{k}.x;
                        reference{k}.y = reference{k}.y;
                        reference{k}.z = reference{k}.z;
                    end
                 end
            else
            end
        else
            % for image model
            if modelDim == 3
                for k = 1:obj.numOfModel
                    if results.movModel
                        oneMPars = obj.exportPars(k,'mPar');          % get mPars
                        img = obj.model{k}.getImage(oneMPars, 'pixelSize',pixelSize,'roiSize',obj.roiSize);

                        % translate lPars to mPars
                        [X,Y,Z] = meshgrid(1:imgSize(1),1:imgSize(2),1:imgSize(3)); 
                        F = griddedInterpolant(img,'cubic', 'nearest');

                        xPos=lPars{k}.x;
                        yPos=lPars{k}.y;
                        zPos=lPars{k}.z;
                        xrot=lPars{k}.xrot;
                        yrot=lPars{k}.yrot;
                        zrot=lPars{k}.zrot;

                        XMean = mean(X(:)); YMean = mean(Y(:)); ZMean = mean(Z(:));
                        cenX = X-XMean;
                        cenY = Y-YMean;
                        cenZ = Z-ZMean;
                        [cenX, cenZ] = rotcoord(cenX(:),cenZ(:),-yrot*pi/180);
                        [cenY, cenZ] = rotcoord(cenY(:),cenZ(:),-xrot*pi/180);
                        [cenX, cenY] = rotcoord(cenX(:),cenY(:),-zrot*pi/180);
                        valNew = F(cenX+XMean-(xPos)./pixelSize, cenY+YMean-(yPos)./pixelSize, cenZ+ZMean-(zPos)./pixelSize);
                        ind = sub2ind(imgSize, X(:),Y(:),Z(:));
                        finalImg(ind)=finalImg(ind)+valNew;
                    else
                        oneMPars = obj.exportPars(k,'mPar');          % get mPars
                        img = obj.model{k}.getImage(oneMPars, 'pixelSize',pixelSize,'roiSize',obj.roiSize);
                        if k == 1
                            finalImg = img;
                            if results.normalizationMax
                                img = img/max(img,1:2);
                            end
                        else
                            oneMPars = obj.exportPars(k,'mPar');          % get mPars
                            oneLPars = obj.exportPars(k,'lPar');          % get mPars
                            img = obj.model{k}.getImage(oneMPars, 'pixelSize',pixelSize,'roiSize',obj.roiSize);
                            if results.normalizationMax
                                img = img/max(img,1:2);
                            end

                            % translate lPars to mPars
                            [X,Y,Z] = meshgrid(1:imgSize(1),1:imgSize(2),1:imgSize(3)); 
                            F = griddedInterpolant(img,'cubic', 'nearest');

                            xPos=oneLPars.x;
                            yPos=oneLPars.y;
                            zPos=-oneLPars.z;
                            xrot=oneLPars.xrot;
                            yrot=oneLPars.yrot;
                            zrot=oneLPars.zrot;

                            XMean = mean(X(:)); YMean = mean(Y(:)); ZMean = mean(Z(:));
                            cenX = X-XMean;
                            cenY = Y-YMean;
                            cenZ = Z-ZMean;
                            [cenX, cenZ] = rotcoord(cenX(:),cenZ(:),-yrot*pi/180);
                            [cenY, cenZ] = rotcoord(cenY(:),cenZ(:),-xrot*pi/180);
                            [cenX, cenY] = rotcoord(cenX(:),cenY(:),-zrot*pi/180);
                            valNew = F(cenX+XMean-(xPos)./pixelSize, cenY+YMean-(yPos)./pixelSize, cenZ+ZMean-(zPos)./pixelSize);
                            ind = sub2ind(imgSize, X(:),Y(:),Z(:));
                            finalImg(ind) = finalImg(ind)+valNew;
                        end
                    end
                end
            else
                for k = 1:obj.numOfModel
                    oneMPars = obj.exportPars(k,'mPar');          % get mPars
                    img = obj.model{k}.getImage(oneMPars, 'pixelSize',pixelSize,'roiSize',obj.roiSize);
                    if results.normalizationMax
                        img = img/max(img,[],1:2);
                    end
                    % translate lPars to mPars
                    [X,Y] = meshgrid(1:imgSize(1),1:imgSize(2)); 

                    F = griddedInterpolant(img,'cubic', 'nearest');

                    xPos=lPars{k}.x;
                    yPos=lPars{k}.y;
                    zrot=lPars{k}.zrot;

                    XMean = mean(X(:)); YMean = mean(Y(:));
                    cenX = X-XMean;
                    cenY = Y-YMean;
                    [cenX, cenY] = rotcoord(cenX(:),cenY(:),-zrot*pi/180);
                    % the X and Y coordinates need not to be inverted using griddedInterpolant
                    valNew = F(cenX+XMean-(xPos)./pixelSize, cenY+YMean-(yPos)./pixelSize);
                    ind = sub2ind(imgSize, X(:),Y(:));
                    finalImg(ind)=finalImg(ind)+valNew;
                end
            end


            switch projection
                case 'xy'
                    finalImg = squeeze(mean(finalImg,3));
                case 'xz'
                    finalImg = squeeze(mean(finalImg,2));
                case 'yz'
                    finalImg = squeeze(mean(finalImg,1));
            end

            if ~isempty(results.axes)
                ax = results.axes;
            else
                fig = figure('Name', 'Fit result');
                ax = axes(fig);
            end
            imagesc(ax, finalImg')
        end
        if isequal(results.plotType, 'point')
                if ~isempty(results.axes)
                    ax = results.axes;
                else
                    fig = figure('Name', 'Fit result');
                    ax = axes(fig);
                end
                for k = 1:obj.numOfModel
                    hold(ax, 'on')
                    plot3(ax, reference{k}.x, reference{k}.y, reference{k}.z, ' ok', 'MarkerSize', 4, 'MarkerFaceColor', 'r')
                end
        end        
        % projection for locs
        if isequal(results.plotType, 'image')
            if ~results.movModel
                locs = obj.locsHandler(locs,lPars{k});
            end
            switch projection
                case 'xy'
                    dim1 = locs.xnm; dim2 = locs.ynm;
                case 'xz'
                    dim1 = locs.xnm; dim2 = locs.znm;
                case 'yz'
                    dim1 = locs.ynm; dim2 = locs.znm;
            end
        end

        if ~isfield(locs,'layer')
            locs.layer = ones(size(locs.xnm));
        end
        layerType = unique(locs.layer);
        col = {'b','g'};

        if isequal(obj.model{whichModel}.modelType, 'discrete')
            hold(ax, 'on')
            for k = length(layerType):-1:1
                lOneLayer = locs.layer == layerType(k);
                if isvarname('lLocsUpRing')
                    plot3(ax, locs.xnm(lLocsUpRing), locs.ynm(lLocsUpRing), locs.znm(lLocsUpRing), ' ow', 'MarkerSize', 8, 'MarkerFaceColor', [0 0 1])
                    plot3(ax, locs.xnm(~lLocsUpRing), locs.ynm(~lLocsUpRing), locs.znm(~lLocsUpRing), ' ow', 'MarkerSize', 8, 'MarkerFaceColor', [0.5 0.5 1])
                else
                    plot3(ax, locs.xnm(lOneLayer), locs.ynm(lOneLayer), locs.znm(lOneLayer), ' ob', 'MarkerSize', 6, 'MarkerFaceColor', col{k})
                end
            end
        else
            hold(ax, 'on')
            for k = length(layerType):-1:1
                lOneLayer = locs.layer == layerType(k);
                if isequal(results.plotType, 'image')
                    plot(ax, (dim1(lOneLayer)+obj.roiSize/2)/pixelSize, (dim2(lOneLayer)+obj.roiSize/2)/pixelSize, ' ow', 'MarkerSize', 3, 'MarkerFaceColor', col{k})
                else
                    plot3(ax, locs.xnm(lOneLayer), locs.ynm(lOneLayer),locs.znm(lOneLayer), ' ob', 'MarkerSize', 6, 'MarkerFaceColor', col{k})
                end
            end
        end
    axis(ax, 'equal')
end