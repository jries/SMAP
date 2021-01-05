classdef cellLevelQC<interfaces.DialogProcessor&interfaces.SEProcessor

    methods
        function obj=cellLevelQC(varargin) 
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={};
        end
        function out=run(obj,p)
            out=runintern(obj,p);
         
        end
        
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function pard=guidef(obj)
mainfile = obj.getPar('mainfile');
if isempty(mainfile)
    defaultDir='';
    defaultFn='';
else
    
 [defaultDir, defaultFn]=fileparts();
end
   % [defaultDir, defaultFn]= fileparts(obj.P.par.mainfile.content); %this creates errors. Somehow this plugin gets loaded in Linux before data is loaded.

    pard.t_shakingCell.object = struct('Style','text','String','Shaking cells (bandwidth&percentile&cutoff(z-score))');
    pard.t_shakingCell.position = [1,1];
    pard.t_shakingCell.Width = 2;

    pard.bandwidth.object = struct('Style','edit','String', 100);
    pard.bandwidth.position = [1,3];
    pard.bandwidth.Width = 0.5;
    
    pard.pretile.object = struct('Style','edit','String',75);
    pard.pretile.position = [1,3.5];
    pard.pretile.Width = 0.5;
    
    pard.zscoreCutoff.object = struct('Style','edit','String',1.5);
    pard.zscoreCutoff.position = [1,4];
    pard.zscoreCutoff.Width = 0.5;
    
    pard.zscoreCutoff2.object = struct('Style','edit','String',2);
    pard.zscoreCutoff2.position = [1,4.5];
    pard.zscoreCutoff2.Width = 0.5;
     
    pard.t_saveTo.object = struct('Style', 'text','String','Save to');
    pard.t_saveTo.position = [2,1];
    pard.t_saveTo.Width = 1;
    
    pard.saveDir.object = struct('Style', 'edit','String',defaultDir);
    pard.saveDir.position = [2,2];
    pard.saveDir.Width = 2;
    
    pard.saveDirSelect.object = struct('Style', 'pushbutton','String','...', 'Callback', {{@callback_folder,obj}});
    pard.saveDirSelect.position = [2,4];
    pard.saveDirSelect.Width = 0.5;
    
    pard.saveFn.object = struct('Style', 'edit','String',defaultFn);
    pard.saveFn.position = [3,2];
    pard.saveFn.Width = 2;

    pard.plugininfo.type='ROI_Analyze';

    pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi','layer1_','layer2_','se_sitepixelsize'};
end

function out=runintern(obj,p)
    out = [];
        %% Init.: getting all the necessary data
    stdClusters = cell([obj.locData.SE.numberOfCells 1]);             % standard deviation of clusters identified using mean-shift clustering
    chNLocs = cell([obj.locData.SE.numberOfCells 1]);                 % number of locs in different channels
    se = obj.locData.SE;
    locs=obj.locData.getloc({'xnm','ynm','clusterdensity','channel','filenumber','locprecnm','locprecnm','PSFxnm'},'position','all','removeFilter',1);
    for k = 1:obj.locData.SE.numberOfCells
        % getting all the necessary data
        fileNumber = obj.locData.SE.cells(k).info.filenumber;
        if k==1||fileNumber ~= fileNumberLast                             % only excuated when file number changed
            indg = locs.filenumber == obj.locData.SE.cells(k).info.filenumber;
            indg=indg&locs.clusterdensity> obj.locData.SE.files(fileNumber).evaluation.densityCutoff;
        end
        %%  Main
        % get all locs in one cell
        
        locsx = locs.xnm(indg);
        locsy = locs.ynm(indg);
        lInCellX = ismember(round(locsx./obj.locData.SE.cells(k).evaluation.preBoundary.pixelSize), obj.locData.SE.cells(k).evaluation.preBoundary.maskPixelSub(:,2));
        lInCellY = ismember(round(locsy./obj.locData.SE.cells(k).evaluation.preBoundary.pixelSize), obj.locData.SE.cells(k).evaluation.preBoundary.maskPixelSub(:,1));
        locsInTheCell = lInCellX&lInCellY;
        locsx = locsx(locsInTheCell);
        locsy = locsy(locsInTheCell);
        
%         % run meanshift clustering and get labels of group for all locs
%         [~,grplabel] = MeanShiftCluster([locsx locsy]',100,0);

%         [~,rankLab] = ismember(grplabel,unique(grplabel));
%         stdAll = splitapply(@(x,y)stdPCAgrp(x,y),locsx,locsy,rankLab');
%         stdAll = grpstats([locsx locsy],grplabel,{'std'});
        
        % summarized std in x and in y 
%         if 0
%             xystd = sqrt(stdAll(:,1).^2.+stdAll(:,2).^2);
%             xystd = xystd(xystd~=0);
%             stdClusters{k} = xystd;
%         else
%             xystd = stdAll(stdAll~=0);
%             stdClusters{k} = xystd;
%         end

        chLabel = locs.channel(indg);
        chLabel = chLabel(locsInTheCell);
        chNLocs{k} = grpstats(chLabel, chLabel, {'numel'});
        fileNumberLast = fileNumber;
    end
    
%     % Get the measurement representing the variation of the clusters in a cell
%     stdClustersPrctile = cellfun(@(x)prctile(x, p.pretile),stdClusters);
%     figure; histogram(stdClustersPrctile,70)
% %     inOrder = stdClustersPrctile90>88.3;
% %     cell2Rm_stdClusters = getFieldAsVector(obj.locData.SE.cells(inOrder), 'ID');
% %     
%     zs = zscore(stdClustersPrctile);
%     figure; histogram(zs, 100);
%     [valStd,rankStd] = sort(zs, 'descend');                                                     % Cells' rank on the zscores
%     inOrder = valStd>p.zscoreCutoff;
%     cell2Rm_stdClusters = getFieldAsVector(obj.locData.SE.cells(rankStd(inOrder)), 'ID');
%     allFrames = cell([length(cell2Rm_stdClusters) 1]);
%     fig = figure(9954);
%     ax = axes(fig);
%     for k = 1:se.numberOfCells
%         allFrames{k,1} = obj.locData.SE.plotcell(obj.locData.SE.cells(rankStd(k)),ax,[]).image;
%     end
%     
%     nrow = ceil(se.numberOfCells/10);
    if 0
        fig = figure;
        a = axes(fig);
        imgSize = obj.getPar('se_cellfov')/obj.getPar('se_cellpixelsize');
        set(fig, 'Position', [0, 0, imgSize, imgSize])
        set(a, 'Position', [0, 0, imgSize, imgSize])
        a.XLim = [0 imgSize];
        a.YLim = [0 imgSize];
        for k = 1:se.numberOfCells
            % add ROIs' ID labels
            text(a, .05,.9,num2str(valStd(k)),'FontSize',38,'FontWeight','bold');
            F = getframe(a,[0, 0, imgSize, imgSize]);
            cla(fig)
            F = F.cdata==0;
            allFrames{k,1}(F==1) = 255;
        end
        close(fig)
        img = montage(allFrames, 'Size', [nrow 10], 'ThumbnailSize', [], 'BackgroundColor', 'white', 'BorderSize', [2 2]);
        imwrite(img.CData, [p.saveDir '\' p.saveFn '_shakingCells.png']);
    end    
   
    % Calutate the ratio of number of locs between the two channels
    chNRatio = cellfun(@(x)x(2)/x(1),chNLocs);
    figure; histogram(chNRatio,100)
    
    if 0
        % perform a guass2 fitting to the kernel density
        [fig,xi] = ksdensity(chNRatio,0:0.01:10, 'Bandwidth',0.05);
        f = fit(xi', fig', 'gauss2');

        % get the curves of the fitting 
        xq = 0:0.05:10;
        val1 = gauss_distribution(xq,f.b1,f.c1);
        val2 = gauss_distribution(xq,f.b2,f.c2);
        figure; plot(xq, f.a1.*val1, xq, f.a2.*val2)

        % detect the intersection of the two guassian distributions in given by the fitting
        [iX,iY] = intersections(xq, f.a1.*val1, xq, f.a2.*val2);

        % detect cells with a chNRatio bellow the threshold determined by the guass2 fitting
        chNRinOrder = chNRatio<iX;
        cell2Rm_chNRatio = getFieldAsVector(obj.locData.SE.cells(chNRinOrder), 'ID');
    else
        zsChNRatio = zscore(chNRatio);
        figure; histogram(zsChNRatio ,100)
        inOrder=zsChNRatio<p.zscoreCutoff2;
        [valChNR,rankChNR] = sort(zsChNRatio);
        cell2Rm_chNRatio = getFieldAsVector(obj.locData.SE.cells(inOrder), 'ID');
        
        figChNR = figure;
        ax = axes(figChNR);
        for k = 1:se.numberOfCells
            allFrames{k,1} = obj.locData.SE.plotcell(obj.locData.SE.cells(rankChNR(k)),ax,[]).image;
        end

        %== 1109: put numbers (measurements) to the figures
        if 0
            fig = figure;
            a = axes(fig);
            set(fig, 'Position', [0, 0, imgSize, imgSize])
            set(a, 'Position', [0, 0, imgSize, imgSize])
            a.XLim = [0 imgSize];
            a.YLim = [0 imgSize];
            for k = 1:se.numberOfCells
                % add ROIs' ID labels
                text(a, .05,.9,num2str(valChNR(k)),'FontSize',38,'FontWeight','bold');
                F = getframe(a,[0, 0, imgSize, imgSize]);
                cla(fig)
                F = F.cdata==0;
                allFrames{k}(F==1) = 255;
            end
            close(fig)
            img = montage(allFrames, 'Size', [nrow 10], 'ThumbnailSize', [], 'BackgroundColor', 'white', 'BorderSize', [2 2]);
            imwrite(img.CData, [p.saveDir '\' p.saveFn '_chNLocsRatio.png']);
        end
    end
    % label the sites with their statuses
%     both = intersect(cell2Rm_stdClusters, cell2Rm_chNRatio);
%     shakingCell = setdiff(cell2Rm_stdClusters, cell2Rm_chNRatio);
%     nonPerCell = setdiff(cell2Rm_chNRatio, cell2Rm_stdClusters);
    nonPerCell = cell2Rm_chNRatio;
    
    sitesParentCells = getFieldAsVector(obj.locData.SE.sites, 'info.cell');
    
%     inShakingCell = ismember(sitesParentCells, shakingCell);
    inNonPerCell = ismember(sitesParentCells, nonPerCell);
%     inBoth = ismember(sitesParentCells, both);
    
%     for k = find(inShakingCell)
%         obj.locData.SE.sites(k).annotation.list3.value=2;
%     end
    for k = find(inNonPerCell)
        obj.locData.SE.sites(k).annotation.list3.value=3;
    end
%     for k = find(inBoth)
%         obj.locData.SE.sites(k).annotation.list3.value=4;
%     end
end

function f = gauss_distribution(x, mu, s)
    p1 = -.5 * ((x - mu)/s) .^ 2;
    p2 = (s * sqrt(2*pi));
    f = exp(p1) ./ p2; 
end

function grpStd = stdPCAgrp(x,y)
    xy = [x y];
    xy = xy - mean(xy,2);
    C = cov(xy);
    [V,D] = eig(C);
    newxy = V * xy';
    newxy = newxy';
    newxy = fliplr(newxy);
    grpStd = std(newxy(:,1));
end

function callback_folder(a,b,obj)
    saveDir = obj.getSingleGuiParameter('saveDir');
    newSaveDir = uigetdir(saveDir);
    obj.setGuiParameters(struct('saveDir', newSaveDir));
end