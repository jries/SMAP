function [ax, v, info] = rotCoordNMkImg(obj, varargin)
% 
% Usage:
%   rotCoordNMkImg(obj, ax, modCoord, locsCoord, rotVizAlt, pixelSize, mode,
% section, lut)
%
% Args:
%   ax:
%   modCoord:
%   locsCoord:
%   rotVizAlt: a vecotr of [az el], where az stands for azimuth and
%   elevation angles.
%   pixelSize:
%   mode:
%   section:
%   lut:

if isa(varargin{1},'matlab.graphics.axis.Axes')
    ax = varargin{1};
    modCoord = varargin{2};
    locsCoord = varargin{3};
    rotVizAlt = varargin{4};
    pixelSize = varargin{5};
    mode = varargin{6};
    section = varargin{7};
    lut = varargin{8};
    varargin = varargin(9:end);
else
    fig = figure;
    ax = axes(fig);
    modCoord = varargin{1};
    locsCoord = varargin{2};
    rotVizAlt = varargin{3};
    pixelSize = varargin{4};
    mode = varargin{5};
    section = varargin{6};
    lut = varargin{7};
    varargin = varargin(8:end);
end

p = inputParser;
p.addParameter('rotMov',false);
p.addParameter('imax',nan);
p.parse(varargin{:});
p = p.Results;

info = [];

% Rotate coordiantes and make an image
forRevY = -1; % this factor was introduced to make the reverse YDir also applies
[locsCoord.xnm,locsCoord.ynm,locsCoord.znm] = rotAzEl(locsCoord.xnm,forRevY.*locsCoord.ynm,locsCoord.znm, rotVizAlt(1), -rotVizAlt(2));
lSectionLocs = locsCoord.znm >= -section & locsCoord.znm <= section;

% Either render an image of the data or the model based on what the user
% selected.
switch mode
    case 'Model'
        v = zeros([obj.roiSize./pixelSize obj.roiSize./pixelSize]);
        % Go through all models (layer-based)
        roiks = 2.7;
        if isempty(obj.getTemp('gausstemplate'))
            G = creategausstemplate(roiks);
            obj.setTemp('gausstemplate',G);
        else
            G = obj.getTemp('gausstemplate');
        end
        
        for k = 1:obj.numOfModel
           
            [modCoord{k}.x,modCoord{k}.y,modCoord{k}.z] = rotAzEl(modCoord{k}.x,modCoord{k}.y,modCoord{k}.z, rotVizAlt(1), -rotVizAlt(2));
            if obj.model{k}.fixSigma
                thisImg = getModelImg(modCoord{k}.x, modCoord{k}.y, 'roiSize', obj.roiSize, 'pixelSize', pixelSize, 'sigma', obj.model{k}.sigma, 'gausstemplate',G,'norm',modCoord{k}.n)';
            else
                variation = obj.getVariable(['par.m' num2str(k) '.lPar.variation']);
                thisImg = getModelImg(modCoord{k}.x, modCoord{k}.y, 'roiSize', obj.roiSize, 'pixelSize', pixelSize, 'sigma', (sqrt(median(locsCoord.locprecnm)^2+variation^2)+obj.model{k}.sigmaFactor(2)).*obj.model{k}.sigmaFactor(1), 'gausstemplate',G,'norm',modCoord{k}.n)';
            end
            thisImg = thisImg/max(thisImg,[],1:2);
            if obj.numOfLayer==1
                thisImg = ind2rgb(round(thisImg.*255./obj.numOfModel), mymakelut('cyan'));
            else
                thisImg = ind2rgb(round(thisImg.*255./obj.numOfModel), mymakelut(obj.model{k}.displayLut));
            end
            v = v + thisImg;
        end
        hold on
        imagesc(ax, v);
        axis(ax,'equal')
        hold off
        
        % Deal with locs
        for k = 1:max(locsCoord.layer)
            lLayer = locsCoord.layer==k;
            if obj.numOfLayer==1
                oneLut = mymakelut(lut{k});
                oneColor = oneLut(140,:);
            elseif k > obj.numOfModel
                oneLut = mymakelut(lut{k});
                oneColor = oneLut(150,:);
            else
                oneLut = mymakelut(obj.model{k}.displayLut);
                oneColor = oneLut(150,:);
            end
            
            hold(ax, 'on');
            vizX = (locsCoord.xnm(lLayer&lSectionLocs)+obj.roiSize/2+pixelSize)./pixelSize;
            vizY = (locsCoord.ynm(lLayer&lSectionLocs)+obj.roiSize/2+pixelSize)./pixelSize;
            plot(ax,  vizX,vizY,' or', 'MarkerEdgeColor','k','MarkerFaceColor',oneColor,'MarkerSize',3.5)
            hold(ax,'off')
        end
        
        % Deal with additional info
        items = obj.getThings2Plot;
        oneItems = items{k};
        for l = 1:length(oneItems)
            oneItem = oneItems(l);
            [oneItem.XData,oneItem.YData,oneItem.ZData] = rotAzEl(oneItem.XData,oneItem.YData,oneItem.ZData, rotVizAlt(1), -rotVizAlt(2));
            oneItem.XData = (oneItem.XData+obj.roiSize/2+pixelSize)./pixelSize;
            oneItem.YData = (oneItem.YData+obj.roiSize/2+pixelSize)./pixelSize;
            oneItem.ZData = (oneItem.ZData+obj.roiSize/2+pixelSize)./pixelSize;
            
            oneItem = rmfield(oneItem,'ZData');
            fn = fieldnames(oneItem);
            hold(ax, 'on');
            h = plot(ax,oneItem.XData,oneItem.YData);
            hold(ax,'off')
            for f = 1:length(fn)
                if ~any(strcmp(fn{f},{'YData','XData'}))
                    set(h,fn{f},oneItem.(fn{f}));
                end
            end
        end
    case 'Data'
        v = zeros([obj.roiSize./pixelSize obj.roiSize./pixelSize]);
        % Deal with locs
        for k = 1:max(locsCoord.layer)
            lLayer = locsCoord.layer==k;
            locsCoordSub = subsetStruct(locsCoord, lLayer&lSectionLocs);
            if isempty(obj.linkedGUI)
                roiks = 2.7;
                if isempty(obj.getTemp('gausstemplate'))
                    G = creategausstemplate(roiks);
                    obj.setTemp('gausstemplate',G);
                else
                    G = obj.getTemp('gausstemplate');
                end
                thisImg = getModelImg(locsCoordSub.xnm, locsCoordSub.ynm, 'roiSize', obj.roiSize, 'pixelSize', pixelSize, 'sigma', locsCoordSub.locprecnm,'gausstemplate',G)';
                % this prevents the pure white image from happening when there
                % is not any locs in the this layer.
                if sum(thisImg(:))>0
                    thisImg = thisImg/max(thisImg,[],1:2);
                end
                maxInt = prctile(thisImg(:),99); % set the saturation
                thisImg(thisImg > maxInt) = maxInt;
                if sum(thisImg(:))>0
                    thisImg = thisImg/max(thisImg,[],1:2);
                end
                thisImg = ind2rgb(round(thisImg.*255./obj.numOfModel), mymakelut(lut{k}));
                v = v + thisImg;
            else
                % use the SMAP render when available
                p_render = obj.linkedGUI.getLayerParameters(k,renderSMAP);
                p_render.sr_pos = [0 0 0];
                p_render.sr_size = repelem(obj.roiSize./2,2);
                p_render.sr_pixrec = pixelSize;
                if p.rotMov
                    p_render.imaxtoggle = 1;
                    imageo = renderSMAP(locsCoordSub, p_render, k);
                    imageo = drawerSMAP(imageo,p_render);
                    imax = imageo.imax*1.8;
                    p_render.imaxtoggle = 0;
                    p_render.imax_min = imax;
                    imageo = renderSMAP(locsCoordSub, p_render, k);
                    imageo = drawerSMAP(imageo,p_render);
                    info.imax = imax;
                else
                    if ~isnan(p.imax)
                        p_render.imax_min = p.imax;
                        p_render.imaxtoggle = 0;
                    end
                    imageo = renderSMAP(locsCoordSub, p_render, k);
                    imageo = drawerSMAP(imageo,p_render);
                end
                thisImg = imageo.image;
                v = v + thisImg;
            end
        end
        imagesc(ax, v);
        axis(ax,'equal')
        
        % Here generate the outline of the model withi the given range of the slice
        if ~isempty(modCoord)
            for k = 1:obj.numOfLayer
                
                [modCoord{k}.x,modCoord{k}.y,modCoord{k}.z] = rotAzEl(modCoord{k}.x,modCoord{k}.y,modCoord{k}.z, rotVizAlt(1), -rotVizAlt(2));
                if rotVizAlt(2)<0
                    lSectionMod = modCoord{k}.z >= -5 & modCoord{k}.z <= 5;
                    modVizX = modCoord{k}.x(lSectionMod);
                    modVizY = modCoord{k}.y(lSectionMod);
                else
                    modVizX = modCoord{k}.x;
                    modVizY = modCoord{k}.y;
                end
                
                %% Get the outline of the model projecion
                % The sorted index idx of the points assiged as the outline is
                % created first for a closed boundary. If the the longest edge
                % are edgeRatioTh-time longer than the second, the edge will be
                % removed to create the open.
                edgeRatioTh = 1.2;
                idx = boundary(double(modVizX), double(modVizY),0);
                if ~isempty(idx)
                    D = sqrt(diff(modVizX(idx)).^2 + diff(modVizY(idx)).^2);
                    [B,I] = sort(D, 'descend');
                    if B(1)/B(2)>edgeRatioTh
                        idx = [idx(I(1)+1:end); idx(1:I(1))];
                        idx = unique(idx,'stable');
                    end
                else
                    idx = 1:length(modVizX);
                end
                
                oneLut = mymakelut(obj.model{k}.displayLut);
                oneColor = oneLut(150,:);
                
                hold(ax, 'on');
                plot(ax,  (modVizX(idx)+obj.roiSize/pixelSize)./2,+(modVizY(idx)+obj.roiSize/pixelSize)./2,'- r', 'Color','w','LineWidth',5)
                hold(ax,'off')
            end
            
        end
end
end
%%
%

function [x,y,z] = rotAzEl(x,y,z, az, el)
% Transform the localization based on the lPars
% Input Arguments
% locs: localizations
% lPars: localization parameters
% need to implement scaling and weight
az = az*pi/180;                      % from degree to radians
el = el*pi/180;
% shift the coordinates
[x,y] = rotcoord(x,y,az);
[y,z] = rotcoord(y,z,el);
end

%%
%

function v = getModelImg(x,y,varargin)
p = inputParser;
p.addParameter('sigma',12)
p.addParameter('pixelSize',5)
p.addParameter('roiSize',500)
p.addParameter('gausstemplate',[])
p.addParameter('norm',[])
p.parse(varargin{:});
results = p.Results;

roiks = 2.7;
sigmax = results.sigma;
roiSize = results.roiSize;
pixelSize = results.pixelSize;
bound = [-roiSize/2 roiSize/2];
rangeSet(1) = ceil(roiSize);
rangeSet(2) = rangeSet(1);
rangeSet(3) = rangeSet(1);
rangeSet = rangeSet./pixelSize;
sig = zeros(size(x));

if isempty(results.norm)
    norm = ones(size(x));
else
    norm = results.norm;
end

G=results.gausstemplate;
v =gaussrender_elliptc((single(x)-bound(1))./pixelSize,(single(y)-bound(1))./pixelSize,uint32(rangeSet(1:2)),single(sig+sigmax/pixelSize),single(sig+sigmax/pixelSize),single(G.template),single(G.sigmatemplate),...
    single(roiks),single(norm),int32(0),single(0),single(0), single([0 0]));
end
%%
%

function gausstemplate=creategausstemplate(roiks) % create template
% global gausstemplate
% sigmatemplate=10;
sizegauss=300;
sigmatemplate=(sizegauss)/(2*roiks)/2; %for 2.5 sigma in both directions
xg=-sizegauss:sizegauss;
[Xg,Yg]=meshgrid(xg,xg);
template=exp(-((Xg).^2+(Yg).^2)/2/sigmatemplate^2);
gausstemplate.template=template;
gausstemplate.sizegauss=sizegauss;
gausstemplate.sigmatemplate=sigmatemplate;
end