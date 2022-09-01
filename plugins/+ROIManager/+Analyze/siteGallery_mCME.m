classdef siteGallery_mCME<interfaces.DialogProcessor&interfaces.SEProcessor
    properties
        fig = gobjects(1);
        fitterGUI_name
    end
    methods
        function obj=siteGallery_mCME(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'se_viewer', 'se_siteroi', 'se_sitefov'};
            obj.showresults=true;
            obj.showloadsavepar=true;
        end
        
        function out=run(obj,p)
            global pan pantext
            pan = [];
            pantext = [];
            % basic info.
            se = obj.locData.SE;
            roiSize = se.P.par.se_siteroi.content;
            sites = se.sites;
            sites2plot = p.sites;
            spBtSites = 1;
            
            % hack an evaluate plug-in in order to use the obj.getLocs(...)
            fdcal=figure(233);
            dcal=plugin('ROIManager','Evaluate','generalStatistics',fdcal,obj.P);
            dcal.attachLocData(obj.SE.locData);
            dcal.makeGui;
            
            % [to-do] allow selection
            fitterGUI_name = 'LocMoFitGUI_2';
            obj.fitterGUI_name = fitterGUI_name;
            %fitterGUI_name = 'LocMoFitGUI';
            eval = obj.locData.SE.processors.eval.guihandles.modules.Data(:,2);
            idxFitterGUI = strcmp(fitterGUI_name,eval);
            fitter = copy(se.processors.eval.processors{idxFitterGUI}.fitter);
            
            obj.setPar('fitter',fitter)
            
            % check the model name (to differentiate mCME and NPC)
            fn = fieldnames(obj.locData.SE.processors.eval.children);
            idx = strcmp(fn, fitterGUI_name);
            modelName = class(obj.locData.SE.processors.eval.processors{idx}.fitter.model{1}.modelObj);
                    
            % parameters
            view = strsplit(p.view,' ');
            layout = p.layout.selection;
            plotModelSchem = p.plotModelSchem;
            dim = p.dim;
            section = p.section;
            
            siteID = getFieldAsVector(se.sites, 'ID');
            if ischar(sites2plot)&&startsWith(sites2plot, '#')
                siteOrder = str2num(sites2plot(2:end));
            else
                [~,siteOrder] = ismember(sites2plot,siteID);
            end
            obj.setPar('siteOrder', siteOrder)
            subSites = se.sites(siteOrder);
            numOfPickedSites = length(subSites);
            obj.setPar('numOfPickedSites',numOfPickedSites);
            
            obj.setPar('dim', dim);
            nSitePage = prod(dim);
                        
            if plotModelSchem
                extraRow = 1;
            else
                extraRow = 0;
            end
            obj.setPar('extraRow', extraRow);
            
            % different views
            numOfView = length(view);
            obj.setPar('numOfView', numOfView);
            
            % layout
            obj.setPar('layout',layout)
            
            switch layout
                case 'vertical'
                    panMode = 'v';
                case 'horizontal'
                    panMode = 'h';
            end
            
            nPage = ceil(numOfPickedSites./nSitePage);
            obj.setPar('nPage',nPage)
            nSiteLastPage = rem(numOfPickedSites,nSitePage);
            if ~isempty(obj.fig)
                delete(obj.fig);
            end
            for np = 1:nPage
                
                obj.fig(np) = figure;
                f = obj.fig(np);
                f.Name = 'siteGallery_mCME_result';
                pan{np} = panel(f);
                pan{np}.pack(dim(2), dim(1));
                pan{np}.de.margin = 0;
                pan{np}.margin = [2 2 2 2];
                
                pantext{np} = panel(f, 'add');
                pantext{np}.pack(dim(2), dim(1));
                pantext{np}.de.margin = 0;
                pantext{np}.margin = [2 2 2 2];
                
                pan{np}.ch.margin = [spBtSites/2 spBtSites/2 spBtSites/2 spBtSites/2];
                ch = pan{np}.ch;
                for j = 1:length(ch)
                    chch = ch{j}.ch;
                    for k = 1:length(chch)
                        chch{k}.margin = [spBtSites/2 spBtSites/2 spBtSites/2 spBtSites/2];
                    end
                end
                
                pantext{np}.ch.margin = [spBtSites/2 spBtSites/2 spBtSites/2 spBtSites/2];
                ch = pantext{np}.ch;
                for j = 1:length(ch)
                    chch = ch{j}.ch;
                    for k = 1:length(chch)
                        chch{k}.margin = [spBtSites/2 spBtSites/2 spBtSites/2 spBtSites/2];
                    end
                end
                % [current]
                f.SizeChangedFcn = {@whenFSizeChanged, f.SizeChangedFcn, {}, @setFSize, {obj}};
                f.CloseRequestFcn = {@closeAllFig,obj};
                
                if np == nPage
                    if nSiteLastPage==0
                        nSiteThisPage = nSitePage;
                    else
                        nSiteThisPage = nSiteLastPage;
                    end
                else
                    nSiteThisPage = nSitePage;
                end
                
                 
                for k = 1:nSiteThisPage
                    if ~strcmp(modelName, 'NPCPointModel_flexible2')
                        plotDirection = 'down';
                    else
                        plotDirection = 'right';
                    end
                    [panR,panC] = tilePosition(k, dim, plotDirection);
                    siteInd = (np-1)*nSitePage+k;
                    panSite = pan{np}(panR,panC);
                    panSite.pack(panMode, numOfView+extraRow);  % for data
                    panSite.de.margin = 0;
                    % Borrow the evaluate plug-in to use getLocs(obj,...)
                    dcal.site=subSites(siteInd);
                    dcal.site.image = se.plotsite(subSites(siteInd));
                    % Check all the possible layers (up to 6)
                    for l = 1:6
                        layercheck = obj.getPar(['layer' num2str(l) '_layercheck']);
                        if layercheck
                            [locsSiteOne,indlocOne] = dcal.getLocs({'xnmrot','ynmrot','znm','locprecnm', 'locprecznm'},'size',roiSize','grouping', 'grouped','layer',l); % per ROI info.
                            if l == 1
                                locsSite = locsSiteOne;
                                indloc = indlocOne;
                                locsSite.layer = ones(size(locsSiteOne.xnm));
                                locsSite.xnm = locsSite.xnmrot;
                                locsSite.ynm = locsSite.ynmrot;
                            else
                                locsSite.layer = [locsSite.layer; ones(size(locsSiteOne.xnm))*l];
                                locsSite.xnm = [locsSite.xnm; locsSiteOne.xnmrot];
                                locsSite.ynm = [locsSite.ynm; locsSiteOne.ynmrot];
                                locsSite.znm = [locsSite.znm; locsSiteOne.znm];
                                locsSite.locprecnm = [locsSite.locprecnm; locsSiteOne.locprecnm];
                                locsSite.locprecznm = [locsSite.locprecznm; locsSiteOne.locprecznm];
                            end
                        end
                    end
                    fitter.allParsArg = subSites(siteInd).evaluation.(fitterGUI_name).allParsArg;
                    fitter.setParArg('m1.lPar.variation', 'value',0);
                    
                    % [to-do] here need to generized so that the models are not
                    % limited to the first one.
                    fitter.model{1}.sigma = p.isoBlurr;
                    fitter.model{1}.fixSigma = true;
                    fitter.model{1}.pixelSize = p.pixelSize;
                    fitter.refPoint_spacing = p.isoGap;
                    fitter.roiSize = 500;
                    [~,modViz] = fitter.plot(locsSite,'plotType','point', 'doNotPlot', true); % get point type visualization
                    lPars = fitter.exportPars(1,'lPar');
                    locsViz = fitter.locsHandler(locsSite, lPars,1);
                    
                    %                 'XYZ' is for model
                    % rotate the view to show the open
                    
                    
                    for indView = 1:numOfView
                        ax = makePlot(fitter, modViz, locsViz, view{indView}, section);
                        set(ax,'YDir','normal')
                        %                     if rot == 1
                        %                        text(ax, page_numOfSites, 30, num2str(subSites(k).ID),'FontSize',50, 'Color','w')
                        %                     end
                        set(ax,'XTick',[], 'YTick', []);
                        set(ax, 'Visible', 'off')
                        axis(ax, 'image')
                        tempFig = ax.Parent;
                        panSite(indView).select(ax);
                        ax.Tag = 'data';
                        close(tempFig);
                    end
                    if plotModelSchem
                        ax = fitter.model{1}.patchPlot(fitter.exportPars(1,'mPar'), 'ele_view', p.tilt, 'pixelSize', p.pixelSize, 'isoCutoff', p.isoCutoff);
                        set(ax,'XLim',[0 fitter.roiSize/p.pixelSize]);
                        set(ax,'XTick',[], 'YTick', [], 'ZTick', []);
                        set(ax, 'XColor', 'none', 'YColor', 'none', 'ZColor', 'none');
                        %                 axis(ax, 'equal')
                        tempFig = ax.Parent;
                        panSite(numOfView+1).select(ax);
                        ax.Tag = 'schematic';
                        close(tempFig);
                    end
                end
                
                %% Export each page
                if p.realTimeSave
                    update_callback([],[],obj,np);
                    drawnow
                    saveOnePage(obj,np);
                end
                
                if np ~= nPage
                    f.Visible = 'off';
                else
                    obj.setPar('currentPage', np)
                end
                
                
            end
            update_callback([],[],obj);
            out = [];
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function ax = makePlot(fitter, modViz, locsViz, view, section)
    if startsWith(view, 'side')
        deg = str2double(view(5:end));
        view = 'side';
    end
    switch view
        case 'top'
            ax = fitter.rotCoordNMkImg(modViz, locsViz, [0 0], 2, 'Data', 500, {'red hot', 'cyan cold'});
        case 'side'
            ax = fitter.rotCoordNMkImg(modViz, locsViz, [deg -90], 2, 'Data', section, {'red hot', 'cyan cold'});
    end
end

function update_callback(a,b,obj, oneNp)
    global pan pantext
    fitterGUI_name = obj.fitterGUI_name;
    se = obj.locData.SE;
    p = obj.getAllParameters;
    fitter = obj.getPar('fitter');
    roiSize = fitter.roiSize;
    pixelSize = p.pixelSize;
    numOfPickedSites = obj.getPar('numOfPickedSites');
    numOfView = obj.getPar('numOfView');
    layout = obj.getPar('layout');
    extraRow = obj.getPar('extraRow');
    dim = obj.getPar('dim');
    siteOrder = obj.getPar('siteOrder');
    subSites = se.sites(siteOrder);
    
    fn = fieldnames(obj.locData.SE.processors.eval.children);
    idx = strcmp(fn, fitterGUI_name);
    modelName = class(obj.locData.SE.processors.eval.processors{idx}.fitter.model{1}.modelObj);
        
    if ~strcmp(modelName, 'NPCPointModel_flexible2')
        plotDirection = 'down';
    else
        plotDirection = 'right';
    end
    
    spBtSites = 1;
    
    figSizeFactor = 2;
    
    if isempty(p.fSize)
        switch layout
            case 'vertical'
                w_fig = figSizeFactor * (dim(2)*(roiSize-p.crop*2)/roiSize+0.1);
                h_fig = figSizeFactor * (dim(1)*(numOfView+extraRow)*(roiSize-p.crop*2)/roiSize+0.1);
            case 'horizontal'
                w_fig = figSizeFactor * (dim(2)*(numOfView+extraRow)*(roiSize-p.crop*2)/roiSize+0.1);
                h_fig = figSizeFactor * (dim(1)*(roiSize-p.crop*2)/roiSize+0.1);
        end
    else
        w_fig = p.fSize(1);
        h_fig = p.fSize(2);
    end
    
    nSitePage = prod(dim);
    nSiteLastPage = rem(numOfPickedSites,nSitePage);
    nPage = obj.getPar('nPage');
    if exist('oneNp','var')
       npRange = oneNp;
    else
       npRange = 1:nPage;
    end
    for np = npRange
        if np == nPage
            if nSiteLastPage==0
                lastID = nSitePage*(np);
            else
                lastID = nSitePage*(np-1)+nSiteLastPage;
            end
        else
            lastID = nSitePage*(np);
        end
                
        
        
        f = pan{np}.figure;
        f.Units = 'centimeters';
        f.Position(3:4) = [w_fig h_fig];
        p1 = findobj(pan{np}.de.axis, 'Tag', 'data');
        p2 = findobj(pan{np}.de.axis, 'Tag', 'schematic');

        set(p1, 'Visible', 'off')
        p1Line = findobj(p1,{'Type','line'},'-and',{'-not',{'Tag','scale bar'}},'-and',{'-not',{'Tag','copy'}});
        p1DashedLine = findobj(p1,'Tag','copy');
        delete(p1DashedLine);
        for ll = 1:length(p1Line)
            oneCopy = copy(p1Line(ll));
            oneCopy.Parent = p1Line(ll).Parent;
            oneCopy.Tag = 'copy';
        end
        p1DashedLine = findobj(p1,'Tag','copy');

        p1SB = findobj(p1,'Tag','scale bar');
        delete(p1SB);
        
    %     axis(p2, 'equal')

        set(p1Line, 'Color', 'k')
        set([p1Line p1DashedLine], 'LineWidth', p.lineWidth)
        set(p1DashedLine, 'Color', p.lineColor, 'LineStyle','--')
        if ~isempty(p2)
            lineObj = findobj(p2, 'type', 'line');
            if isempty(lineObj)
                p2 = reshape(p2, dim(2), [])';
                p2 = p2(:);
                for k = 1:length(p2)
                    [panR,panC] = ind2sub(dim,k);
                    outline_Mod = copy(findobj(pan{np}(panR,panC,numOfView).axis,'Color',p.lineColor));
                    outline_Mod.ZData = outline_Mod.YData;
                    outline_Mod.YData = zeros(size(outline_Mod.YData))+50;
                    outline_Mod.Parent = p2(k);
                end
                axis(p2, 'normal')
                axis(p2, 'equal')
                axis(p2, 'vis3d')
            else
                delete(lineObj)
                axis(p2, 'equal')
                axis(p2, 'vis3d')
            end
        end
        axis(p1, 'image')

    %     axis(p2, 'equal')
        set(p1, 'XLim', [p.crop/pixelSize (roiSize-p.crop)/pixelSize]);
        set(p1, 'YLim', [p.crop/pixelSize (roiSize-p.crop)/pixelSize]);
        set(p2, 'XLim', [p.crop/pixelSize (roiSize-p.crop)/pixelSize]);
    %     axis(p2, 'tight')
        ax = pan{np}(1,1,1).axis;
        addScalebar(ax,'bottom-left', [40 40]./p.pixelSize,100/p.pixelSize);
% 
%         pan{np}.ch.margin = [spBtSites/2 spBtSites/2 spBtSites/2 spBtSites/2];
%         ch = pan{np}.ch;
%         for j = 1:length(ch)
%             chch = ch{j}.ch;
%             for k = 1:length(chch)
%                 chch{k}.margin = [spBtSites/2 spBtSites/2 spBtSites/2 spBtSites/2];
%             end
%         end
%         findObj()
%         labelOn
        
        if ~strcmp(modelName, 'NPCPointModel_flexible2')
            firstLastID = [nSitePage*(np-1)+1 lastID];
            if any(p.labelOrder>0)
                for s = firstLastID(1):firstLastID(2)
                    labels = {'ID','theta','curvature','radius','area'};
                    usedLabels = labels(p.labelOrder>0);
                    usedLabels((p.labelOrder(p.labelOrder>0))) = usedLabels;
                    siteLabel = getLabels(subSites, usedLabels, s, fitterGUI_name);
                    posInPage = rem(s,nSitePage);
                    if posInPage==0
                        posInPage = nSitePage;
                    end
                    [panR,panC] = tilePosition(posInPage, dim, plotDirection);
                    onePantext = pantext{np}(panR,panC);
                    delete(onePantext.axis)
                    onePantext.select(axes('Color','none'))
                    text(onePantext.axis,0.3, 0.05,siteLabel,'Color',[1 1 1],'VerticalAlignment','baseline')
                    hold on
                    text(onePantext.axis,0.01,0.85,num2str(siteOrder(s),'%04.f'),'Color',[1 1 1],'VerticalAlignment','baseline')
                    hold off
                end
            end
        end
        axis(pantext{np}.de.axis,'off');
        hAllText = findobj(pantext{np}.de.axis,'type','text');
        if p.labelOn
            textVisible = 'on';
        else
            textVisible = 'off';
        end
        set(hAllText,'Visible',textVisible)
    end
end

function siteLabel = getLabels(subSites, labels, siteInd, fitterGUI_name)
numOfUsedLabels = length(labels);
for k = 1:numOfUsedLabels
    switch labels{k}
        case 'ID'
            oneLabel = ['Site ' num2str(siteInd)];
        case 'theta'
            val = subSites(siteInd).evaluation.(fitterGUI_name).fitInfo.derivedPars{1}.closingAngle_pub;
            oneLabel = ['\theta = '  num2str(val, '%.1f') char(176)];
        case 'curvature'
            val = subSites(siteInd).evaluation.(fitterGUI_name).fitInfo.derivedPars{1}.curvature;
            oneLabel = ['\it1/R\rm = '  num2str(val, '%.1e') ' nm^{-1}'];
        case 'radius'
            val = subSites(siteInd).evaluation.(fitterGUI_name).fitInfo.derivedPars{1}.radius;
            oneLabel = ['\itR\rm = '  num2str(val, '%.0f') ' nm'];
        case 'area'
            val = subSites(siteInd).evaluation.(fitterGUI_name).fitInfo.derivedPars{1}.realSurfaceArea;
            oneLabel = ['\itA\rm = '  num2str(val, '%.0f') ' nm^2'];
    end
    if k == 1
        siteLabel = oneLabel;
    else
        siteLabel = [siteLabel ', ' oneLabel];
    end
end
end

function whenFSizeChanged(a,b,varargin)
    for k = 1:length(varargin)/2
        varargin{(2*k-1)}(a,b,varargin{2*k}{:});
    end
end

function setFSize(a,b,obj)
    fSize = a.Position(3:4);
	setGuiParameters(obj, struct('fSize',num2str(fSize)))    
end

function pard=guidef(obj)
col2 = 1.8;
col4 = 3.8;

rowRun = 1;
pard.t1.object=struct('String','Site IDs','Style','text');
pard.t1.position=[rowRun,1];
pard.t1.Width=1;

pard.sites.object=struct('String','','Style','edit');
pard.sites.position=[rowRun,col2];
pard.sites.Width=1;

pard.t_pixelSize.object=struct('String','Pixel size','Style','text');
pard.t_pixelSize.position=[rowRun+1,1];
pard.t_pixelSize.Width=1;

pard.pixelSize.object=struct('String','2','Style','edit');
pard.pixelSize.position=[rowRun+1,col2];
pard.pixelSize.Width=1;

pard.t_isoBlurr.object=struct('String','isoBlurr','Style','text');
pard.t_isoBlurr.position=[rowRun+2,1];
pard.t_isoBlurr.Width=1;

pard.isoBlurr.object=struct('String','3','Style','edit');
pard.isoBlurr.position=[rowRun+2,col2];
pard.isoBlurr.Width=1;
pard.isoBlurr.Tooltip = 'The gaussian sigma for blurring. This controls the smoothness of the isosurface model rendering.';

pard.t_isoGap.object=struct('String','isoGap','Style','text');
pard.t_isoGap.position=[rowRun+3,1];
pard.t_isoGap.Width=1;

pard.isoGap.object=struct('String','0.14','Style','edit');
pard.isoGap.position=[rowRun+3,col2];
pard.isoGap.Width=1;
pard.isoGap.Tooltip = 'Gap between sampled points. This controls the sampling rate of the isosurface model rendering.';

pard.t_isoCutoff.object=struct('String','isoCutoff','Style','text');
pard.t_isoCutoff.position=[rowRun+4,1];
pard.t_isoCutoff.Width=1;

pard.isoCutoff.object=struct('String','4','Style','edit');
pard.isoCutoff.position=[rowRun+4,col2];
pard.isoCutoff.Width=1;
pard.isoCutoff.Tooltip = 'Intensity cutoff of the original image for the isosurface model rendering.';

pard.t_tilt.object=struct('String','Tilt','Style','text');
pard.t_tilt.position=[rowRun+5,1];
pard.t_tilt.Width=1;

pard.tilt.object=struct('String','0','Style','edit');
pard.tilt.position=[rowRun+5,col2];
pard.tilt.Width=1;
pard.tilt.Tooltip = 'Tilt angle of the 3d rendering';

pard.set1.object=struct('String','Use','Style','pushbutton','Callback',{{@set1_callback,obj}});
pard.set1.position=[rowRun,3];
pard.set1.Width=0.5;
pard.set1.Tooltip = 'Pre-defined set 1';

pard.set2.object=struct('String','Vesicle','Style','pushbutton','Callback',{{@set2_callback,obj}});
pard.set2.position=[rowRun,3.5];
pard.set2.Width=0.5;
pard.set2.Tooltip = 'Pre-defined set 2';

pard.set3.object=struct('String','Set3','Style','pushbutton','Callback',{{@set3_callback,obj}});
pard.set3.position=[rowRun,4];
pard.set3.Width=0.5;
pard.set3.Tooltip = 'Pre-defined set 3';

pard.t_view.object=struct('String','View','Style','text');
pard.t_view.position=[rowRun+1,3];
pard.t_view.Width=1;

pard.view.object=struct('String','top side0 side45 side90 side135','Style','edit');
pard.view.position=[rowRun+1,col4];
pard.view.Width=1;

pard.t_layout.object=struct('String','Layout','Style','text');
pard.t_layout.position=[rowRun+2,3];
pard.t_layout.Width=1;

pard.layout.object=struct('String','horizontal|vertical','Value',1,'Style','popupmenu');
pard.layout.position=[rowRun+2,col4];
pard.layout.Width=1;

pard.plotModelSchem.object=struct('String','Model schematic','Value',0,'Style','checkbox');
pard.plotModelSchem.position=[rowRun+3,3];
pard.plotModelSchem.Width=1;

pard.t_dim.object=struct('String','Dimension (site)','Style','text');
pard.t_dim.position=[rowRun+4,3];
pard.t_dim.Width=1;

pard.dim.object=struct('String','14 2','Style','edit');
pard.dim.position=[rowRun+4,col4];
pard.dim.Width=1;

pard.t_section.object=struct('String','Section (sideview)','Style','text');
pard.t_section.position=[rowRun+5,3];
pard.t_section.Width=1;

pard.section.object=struct('String',30,'Style','edit');
pard.section.position=[rowRun+5,col4];
pard.section.Width=1;

rowUpdate = 8;

pard.update.object=struct('String','update','Style','pushbutton','callback',{{@update_callback,obj}});
pard.update.position=[rowUpdate+0.5,3.5];
pard.update.Width=1;

pard.up.object=struct('String','<','Style','pushbutton','callback',{{@turnPage_callback,obj}});
pard.up.position=[rowUpdate+1.5,3.5];
pard.up.Width=0.5;

pard.down.object=struct('String','>','Style','pushbutton','callback',{{@turnPage_callback,obj}});
pard.down.position=[rowUpdate+1.5,4];
pard.down.Width=0.5;

pard.saveAll.object=struct('String','Save all','Style','pushbutton','callback',{{@saveFigs,obj}});
pard.saveAll.position=[rowUpdate+2.5,3.5];
pard.saveAll.Width=1;

pard.realTimeSave.object=struct('String','','Style','checkbox','value',1);
pard.realTimeSave.position=[rowUpdate+3.5,3.5];
pard.realTimeSave.Width=0.3;
pard.realTimeSave.Tooltip='Check this box if you want to save a page once it is done.';

pard.saveTo.object=struct('String','Save to','Style','pushbutton','callback',{{@saveTo_callback,obj}});
pard.saveTo.position=[rowUpdate+3.5,3.7];
pard.saveTo.Width=0.5;

pard.saveToPath.object=struct('String','','Style','edit');
pard.saveToPath.position=[rowUpdate+3.5,4.2];
pard.saveToPath.Width=0.8;

pard.t_crop.object=struct('String','Crop','Style','text');
pard.t_crop.position=[rowUpdate,1];
pard.t_crop.Width=1;

pard.crop.object=struct('String','90','Style','edit');
pard.crop.position=[rowUpdate,col2];
pard.crop.Width=1;

pard.t_lineColor.object=struct('String','Line color','Style','text');
pard.t_lineColor.position=[rowUpdate+1,1];
pard.t_lineColor.Width=1;

pard.lineColor.object=struct('String','#0247ff','Style','edit');
pard.lineColor.position=[rowUpdate+1,col2];
pard.lineColor.Width=1;

pard.t_lineWidth.object=struct('String','Line width','Style','text');
pard.t_lineWidth.position=[rowUpdate+2,1];
pard.t_lineWidth.Width=1;

pard.lineWidth.object=struct('String','1.5','Style','edit');
pard.lineWidth.position=[rowUpdate+2,col2];
pard.lineWidth.Width=1;

pard.t_fSize.object=struct('String','Figure size','Style','text');
pard.t_fSize.position=[rowUpdate+3,1];
pard.t_fSize.Width=1;

pard.t_label.object=struct('String','ID/theta/curvature/radius/area','Style','text');
pard.t_label.position=[rowUpdate+4,1];
pard.t_label.Width=1.5;
pard.t_label.Tooltip='The order of the types of labels. 0 for disable. "1 2 0 0 3" for ID/theta/area.';

pard.labelOrder.object=struct('String','1 2 0 0 3','Style','edit');
pard.labelOrder.position=[rowUpdate+4,2.3];
pard.labelOrder.Width=0.8;
pard.labelOrder.Tooltip='The order of the types of labels. 0 for disable. "1 2 0 0 3" for ID/theta/area.';

pard.labelOn.object=struct('String','','Style','checkbox','value',1);
pard.labelOn.position=[rowUpdate+4,3.1];
pard.labelOn.Width=0.3;

pard.fSize.object=struct('String','','Style','edit');
pard.fSize.position=[rowUpdate+3,col2];
pard.fSize.Width=1;

% pard.t2.object=struct('String','Which fitter','Style','text');
% pard.t2.position=[3,1];
% pard.t2.Width=1;
% 
% eval = obj.locData.SE.processors.eval.guihandles.modules.Data(:,2);
% lFitter = startsWith(eval,'LocMoFitGUI');
% options = eval(lFitter);
% 
% pard.fitter.object=struct('String',{options},'value',1,'Style','popupmenu');
% pard.fitter.position=[3,2];
% pard.fitter.Width=1;

pard.plugininfo.type='ROI_Analyze';

end

function set1_callback(a,b,obj)
    se = obj.locData.SE;
    sites = se.sites;
    ID = getFieldAsVector(sites,'ID');
    use = getFieldAsVector(sites,'annotation.use');
    setGuiParameters(obj, struct('sites',num2str(ID(use))))
end

function set2_callback(a,b,obj)
    se = obj.locData.SE;
    sites = se.sites;
    realCloseAngle = [];
    for k = se.numberOfSites:-1:1
        realCloseAngle(k) = sites(k).evaluation.(fitterGUI_name).fitInfo.derivedPars{1}.realCloseAngle+90;
    end
    ID = getFieldAsVector(sites,'ID');
    use = getFieldAsVector(sites,'annotation.use');
    vesicles = ID(realCloseAngle>135&use);
    setGuiParameters(obj, struct('sites',num2str(vesicles)))
end

function set3_callback(a,b,obj)
    disp('Set3 yet to be defined.')
end

function turnPage_callback(a,b,obj)
    global pan
    currentPage = obj.getPar('currentPage');
    if strcmp(a.String, '>') 
        gap = 1;
    elseif strcmp(a.String, '<')
        gap = -1;
    end
    fig_old = pan{currentPage}.figure;
    fig_new = pan{currentPage+gap}.figure;

    fig_old.Visible = 'off';
    fig_new.Visible = 'on';
    fig_new.Position = fig_old.Position;
    obj.setPar('currentPage', currentPage+gap)
end

function closeAllFig(a,b,obj)
    global pan
    for k = 1:length(pan)
        delete(pan{k}.figure)
    end
    delete(a)
end

function saveFigs(a,b,obj)
    [file,path] = uiputfile({'*.pdf;*.png'},'Save all figures to...');
    [~,fName,fExt] = fileparts(file);
    global pan
    update_callback([],[],obj);
    for k = 1:length(pan)
        exportgraphics(pan{k}.figure,[path fName '_' num2str(k) fExt],'ContentType','vector', 'Resolution', 600)
    end
end

function saveOnePage(obj,k)
    global pan
    p = obj.getGuiParameters;
    [path,fName,fExt] = fileparts(p.saveToPath);
    if isempty(fExt)
        fExt = '.png';
    end
    exportgraphics(pan{k}.figure,[path filesep fName '_' num2str(k) fExt],'ContentType','vector', 'Resolution', 600)
end

function saveTo_callback(a,b,obj)
    [file,path] = uiputfile({'*.pdf;*.png'},'Save all figures to...');
    [~,fName,fExt] = fileparts(file);
    p.saveToPath = [path fName];
    obj.setGuiParameters(p)
end