classdef LocMoFit_manager < handle
    % :class:`LocMoFit_manager` a mamager of LocMoFit objects.
    %
    % Prop:
    %   plotSettings: a struct with plot names defined as fieldnames. Each
    %   files is a cell with two elements of variable IDs.
    %   filteringRule: 
    properties
        data                % an SEsites array
        currentData = 1;
        lUseOnly = true;
        usedSites           % define which sites to include here
        grpRule
        plotSettings = [];
        filteringRule
        variableTableCol = [];
        alignSettings       % Converter for alignment.
        handles             % Place to save graphic handles
    end
    properties (Transient)
        parentObj
        sites_locs
        idxCurrentSite         % for registration
    end
    properties (Dependent)
        numOfSites
        numOfFitter
        fitterNames
        useSites
        IDSites
        fileNumberSites
        
        filtering
        
        grp
        
        numOfCh
        NLayer2
        
        variableTable
        convertTable;
    end
    methods
        function obj = LocMoFit_manager(sites, fitter)
            obj.data(1).sites = sites;
            obj.data(1).fitter = fitter;
        end
        function [val, ID] = getVariable(obj, fitterAndID)
            % Only get varibles of use sites
            [val, ID] = getVariable_allSites(obj, fitterAndID);
            if ~isempty(obj.usedSites)
                val = val(obj.usedSites);
            elseif obj.lUseOnly
                val = val(obj.useSites);
            end
        end
        function [val, ID] = getVariable_allSites(obj, fitterAndID)
            IDParts = split(fitterAndID,'.');
            mainFitterID = IDParts{1};
            nameOrID = [IDParts{2} '.' IDParts{3}];
            for k = obj.numOfSites:-1:1
                if isfield(obj.data(obj.currentData).sites(k).evaluation.(mainFitterID), 'allParsArg')
                    fitter = obj.loadInfo(mainFitterID, k);
                    oneVal = fitter.getVariable(nameOrID);
                    if isempty(oneVal)
                        val(k)=nan;
                        ID = [];
                    else
                        [val(k),ID] = fitter.getVariable(nameOrID);
                    end
                else
                    val(k)=nan;
                    ID = [];
                end
                
            end
        end
        function IDs = variableIDs(obj, varargin)
            %% 
            IDs = {};
            for k = 1:obj.numOfFitter
                fitterName = obj.fitterNames{k};
                IDs = [IDs; strcat([fitterName '.'], obj.data(obj.currentData).fitter.(fitterName).getAllParId([],'type','all'))];
            end
        end
        
        function hList = plot(obj)
            % based on settings
            % histogram for 1D data
            % integrated scatter plot for 2D data
            fn = fieldnames(obj.plotSettings);
            lKept = obj.filtering;
            for k = 1:length(fn)
                onePlot = obj.plotSettings.(fn{k});
                if length(onePlot) == 2
                    xData = obj.getVariable(onePlot{1});
                    yData = obj.getVariable(onePlot{2});
                    xData = xData(lKept);
                    yData = yData(lKept);
                    grp = obj.grp(obj.useSites);
                    grp = grp(lKept);
                    fig = figure;
                    t = tiledlayout(fig, 3,3);
                    ax = nexttile(1,[2 2]);
                    ax_2 = nexttile(7,[1 2]);
                    ax_3 = nexttile(3,[2 1]);
                    
                    allPoints = [xData; yData*2000]';
                    
                    grpId = unique(grp);
                    density = zeros(size(grp));
                    for l = 1:length(grpId)
                        density(grp==grpId(l)) = ...
                            countneighbours_versatile(allPoints(grp==grpId(l),:),...
                            allPoints(grp==grpId(l),:),5,0);
                    end
                    if isempty(obj.parentObj)
                        grpScatter(ax, xData, yData, grp,3,-density, 'filled');
                    else
                        grpScatter(ax, xData, yData, grp,3,-density, 'filled');
%                         plotSElink(ax, xData, yData, obj.IDSites(obj.useSites), obj.parentObj.SE,' ko');
                    end
                    
                    % x-axis
                    histogram(ax_2, xData);
                    ax_2.YDir = 'reverse';
                    xlabel(ax_2, setIDFormat(onePlot{1}, 'meaning'))
                    % y-axis
                    histogram(ax_3, yData);
                    
                    % assigning identifiers
                    ax.Tag = 'scatter';
                    ax_2.Tag = 'xHist';
                    ax_3.Tag = 'yHist';
%                     ax_3.YDir = 'reverse';
                    ax_3.View = [90 -90];
                    ax_2.XLim = ax.XLim;
                    ax_3.XLim = ax.YLim;
                    xlabel(ax_3, setIDFormat(onePlot{2}, 'meaning'))
                    ax_3.XAxisLocation = 'top';
                    
%                     tabName = [setIDFormat(onePlot{1}, 'short') '__' setIDFormat(onePlot{2}, 'short')];
                    hList.(fn{k}) = t;
                elseif length(onePlot) == 1
                    data = obj.getVariable(onePlot{1});
                    data = data(lKept);
                    fig = figure;
                    hList.(fn{k}) = axes(fig);
                    histogram(hList.(fn{k}),data,20);
                    title(onePlot{1});
                end
            end
        end
        
        function lFilter = filter(obj, varargin)
            % Use pair parameters like filter('parID',[0 inf])
            % It gets filter for all sites.
            parIDs = varargin(1:2:end);
            parFilter = varargin(2:2:end);
            for k = 1:length(parIDs)
                % Don't change this to obj.getVariable. This function is
                % designed to applied to all sites including the ones
                % marked as not used in the roi manager.
                val = obj.getVariable_allSites(parIDs{k});
                if k == 1
                    lFilter = ones(size(val));
                end
                lOneFilter = val>=parFilter{k}(1);
                lOneFilter = lOneFilter&val<=parFilter{k}(2);
            end
            lFilter = lFilter&lOneFilter;
        end
        
        function addFiltering(obj, varargin)
            % Use pair parameters like obj.addFilter('parID',[0 inf])
            if isempty(obj.filteringRule)
                obj.filteringRule = {};
            end
            parIDs = varargin(1:2:end);
            parFilter = varargin(2:2:end);
            if ~isempty(obj.filteringRule)
                [indExist, locInList] = ismember(parIDs, obj.filteringRule(:,1));
            else
                indExist = 0;
                locInList = 0;
            end
            
            if any(indExist)
                obj.filteringRule(locInList,:) = [parIDs(indExist) parFilter(indExist)];
            end
            
            if any(~indExist)
                n = sum(~indExist);
                obj.filteringRule(end+1:end+n,:) = [parIDs(~indExist) parFilter(~indExist)];
            end
%             obj.filteringRule{} = parIDs
        end
        
        function addGrpRule(obj, varargin)
            % Usage:
            %   obj.addFilter(grp, rule)            
            if isempty(obj.grpRule)
                obj.grpRule = {};
            end
            allGrp = varargin(1:2:end);
            allCond = varargin(2:2:end);
            if ~isempty(obj.grpRule)
                [indExist, locInList] = ismember(allGrp, obj.grpRule(:,1));
            else
                indExist = 0;
                locInList = 0;
            end
            
            if any(indExist)
                obj.grpRule(locInList,:) = [allGrp(indExist) allCond(indExist)];
            end
            
            if any(~indExist)
                n = sum(~indExist);
                obj.grpRule(end+1:end+n,:) = [allGrp(~indExist) allCond(~indExist)];
            end
%             obj.filteringRule{} = parIDs
        end
                
        function batchFit(obj)
        end
        
        function numOfSites = get.numOfSites(obj)
            numOfSites = length(obj.data(obj.currentData).sites);
        end
        
        function numOfFitter = get.numOfFitter(obj)
            numOfFitter = length(fieldnames(obj.data(obj.currentData).fitter));
        end
        
        function fitterNames = get.fitterNames(obj)
            fitterNames = fieldnames(obj.data(obj.currentData).fitter);
        end
        
        function l = get.useSites(obj) 
            l = getFieldAsVector(obj.data(obj.currentData).sites, 'annotation.use');
        end
        
        function ID = get.IDSites(obj) 
            ID = getFieldAsVector(obj.data(obj.currentData).sites, 'ID');
        end
        
        function ID = get.fileNumberSites(obj) 
            ID = getFieldAsVector(obj.data(obj.currentData).sites, 'info.filenumber');
        end
        
        function numOfCh = get.numOfCh(obj)
            numOfCh = length(obj.data(obj.currentData).sites(1).evaluation.generalStatistics.Nlayers);
        end
        function NLayer2 = get.NLayer2(obj) 
            NLayer2 = getFieldAsVectorInd(obj.data(obj.currentData).sites, 'evaluation.generalStatistics.Nlayers', 2);
        end
        
        function convertTable = get.convertTable(obj)
            if ~isempty(obj.parentObj)&&~isempty(obj.parentObj.getPar(''))
            else
            end
        end
        
        function lFilter = get.filtering(obj)
            lFilter = obj.filter(obj.filteringRule{:});
        end
        
        function grp = get.grp(obj)
            grp = zeros(obj.numOfSites,1);
            if ~isempty(obj.grpRule)
                grpInd = [obj.grpRule{:,1}];
                for k = 1:size(grpInd,1)
                    grpId = obj.grpRule{k,1};
                    rule = obj.grpRule{k,2};
                    ruleExprs = regexprep(rule, '(LocMoFitGUI(_\d+)?\.m\d+\.\w+)', ['obj\.' char("getVariable_allSites(\'$1\')")]);          % Parameters
                    grp(eval(ruleExprs)) = grpId;
                end
            end
%             grp = 
        end
               
        function variableTable = get.variableTable(obj)
            usedSites = logical(obj.useSites);
            variableTable = obj.IDSites(usedSites)';
            variableTable(:,2) = obj.fileNumberSites(usedSites)';
            if obj.numOfCh==2
                variableTable(:,3) = obj.NLayer2(usedSites)';
                preColName = ["ID" "filenumber" "NLayer2"];
            else
                preColName = ["ID" "filenumber"];
            end
            numPreCol = size(variableTable,2);
            for k = 1:length(obj.variableTableCol)
                oneID = obj.variableTableCol{k};
                variableTable(:,k+numPreCol) = obj.getVariable(oneID);
            end
            
            % add col names
            colHead = setIDFormat(obj.variableTableCol, 'meaning');
            variableTable = string(variableTable);
            variableTable = [[preColName string(colHead)]; variableTable];
        end
               
        function fitter = loadInfo(obj, mainFitterID, k)
            fitter = obj.data(obj.currentData).fitter.(mainFitterID);
            fitter.fitInfo = obj.data(obj.currentData).sites(k).evaluation.(mainFitterID).fitInfo;
            if isfield(obj.data(obj.currentData).sites(k).evaluation.(mainFitterID), 'externalInfo')
                fitter.externalInfo = obj.data(obj.currentData).sites(k).evaluation.(mainFitterID).externalInfo;
            else
                fitter.externalInfo = 0;
            end
            fitter.allParsArg = obj.data(obj.currentData).sites(k).evaluation.(mainFitterID).allParsArg;
            fitter.getDerivedPars;
        end
        
        function h = createPlot(obj, dataPlot)
            % plot histogram and/or scatter plot based on the settings
            fn = fieldnames(dataPlot);
            for k = 1:length(fn)
                if length(dataPlot.(fn{k}))==2
                    h.(fn{k}) = histogram(dataPlot.(fn{k}));
                    h.(fn{k}) = histogram(dataPlot.(fn{k}));
                else
                    h.(fn{k}) = histogram(dataPlot.(fn{k}));
                end
            end
        end
        
        function exportLocs = registerSites(obj, varargin)
            % to-do: LocMoFitGUI_2 here has to be replaced with
            % obj.currentFitter
            
            p = inputParser;
            fitter = obj.data(obj.currentData).fitter.LocMoFitGUI_2;
            roiSize = obj.parentObj.getPar('se_siteroi');
            p.addParameter('spatialOffset', 0);
            p.addParameter('spatialTrimXY', zeros(1,2));
            p.addParameter('distFromOrigin', 1000);
            p.addParameter('firstBin', false);
            p.parse(varargin{:});
            results = p.Results;
            
            parentObj = obj.parentObj;
            %% alignment           
            if isempty(obj.idxCurrentSite)
                error('No ID of the current site specified.')
                return
            end
            se = obj.parentObj.SE;
            sites = obj.parentObj.SE.sites;
            
            % hack an evaluate plug-in in order to use the obj.getLocs(...)
            fdcal = obj.parentObj.getPar('fdcalForManager');
            if isempty(fdcal)||~isvalid(fdcal)
                fdcal=figure(233);
                obj.parentObj.setPar('fdcalForManager', fdcal)
            end

            dcal=plugin('ROIManager','Evaluate','generalStatistics',fdcal,obj.parentObj.P);
            dcal.attachLocData(se.locData);
            dcal.makeGui;
            
            % Borrow the evaluate plugin to use getLocs(obj,...)
            dcal.site=sites(obj.idxCurrentSite);
            dcal.site.image = se.plotsite(dcal.site.ID);
            [newlocs,indNewLocs] = dcal.getLocs({'xnmrot','ynmrot','znm','locprecnm', 'locprecznm'},'size',roiSize','grouping', 'ungrouped'); % per ROI info.
            
            if results.firstBin
                settings = obj.firstBin_settings;
                % reset the rules
                oldRules = fitter.converterRules;
                fitter.converterRules = [];
                fitter.converterUserDefined = [];

                for k = 1:length(settings.converter)
                    fitter.converter([], settings.converter(k).rule, settings.converter(k).target)
                end
            end
            
            % convert
            % load allParsArg of the specific site
            fitter.allParsArg = sites(obj.idxCurrentSite).evaluation.LocMoFitGUI_2.allParsArg;
            fitter.fitInfo = sites(obj.idxCurrentSite).evaluation.LocMoFitGUI_2.fitInfo;
            fitter.externalInfo = sites(obj.idxCurrentSite).evaluation.LocMoFitGUI_2.externalInfo;
            fitter.convertNow
            
            newlocs.xnm = newlocs.xnmrot;
            newlocs.ynm = newlocs.ynmrot;
            newlocs = fitter.locsRegister(newlocs, fitter.exportPars(1,'lPar'), 1);
            
            if results.firstBin
                fitter.converterRules = oldRules;
            end
            
            % save to the fields
            idxNewLocs = find(indNewLocs);
            
            % Trim locs
            indKept = abs(newlocs.xnm)<(roiSize/2 - results.spatialTrimXY(1));
            indKept = indKept&abs(newlocs.ynm)<(roiSize/2 - results.spatialTrimXY(2));
            fn = fieldnames(newlocs);
            for k = 1:length(fn)
                newlocs.(fn{k}) = newlocs.(fn{k})(indKept);
            end
            
            if results.firstBin
                parentObj.locData.loc.(['xnmaligned_' 'masterAvgMod'])(idxNewLocs) = 0;
                parentObj.locData.loc.(['ynmaligned_' 'masterAvgMod'])(idxNewLocs) = 0;
                parentObj.locData.loc.(['znmaligned_' 'masterAvgMod'])(idxNewLocs) = 0;
            end
            
            idxNewLocs = idxNewLocs(indKept);
            if nargout>0&&results.firstBin
                exportLocs.xnm = newlocs.xnm+results.spatialOffset+results.distFromOrigin;
                exportLocs.ynm = newlocs.ynm;
                exportLocs.znm = newlocs.znm;
            else
                if results.firstBin
                    parentObj.locData.loc.(['xnmaligned_' 'masterAvgMod'])(idxNewLocs) = newlocs.xnm+results.spatialOffset+results.distFromOrigin;
                    parentObj.locData.loc.(['ynmaligned_' 'masterAvgMod'])(idxNewLocs) = newlocs.ynm;
                    if isfield(newlocs,'znm')
                        parentObj.locData.loc.(['znmaligned_' 'masterAvgMod'])(idxNewLocs) = newlocs.znm;
                    end
                else
                    lWithoutVal = parentObj.locData.loc.(['xnmaligned_' 'masterAvg'])(idxNewLocs) == 0;
                    parentObj.locData.loc.(['xnmaligned_' 'masterAvg'])(idxNewLocs(lWithoutVal)) = newlocs.xnm(lWithoutVal)+results.spatialOffset+results.distFromOrigin;
                    parentObj.locData.loc.(['ynmaligned_' 'masterAvg'])(idxNewLocs(lWithoutVal)) = newlocs.ynm(lWithoutVal);
                    if isfield(newlocs,'znm')
                        parentObj.locData.loc.(['znmaligned_' 'masterAvg'])(idxNewLocs(lWithoutVal)) = newlocs.znm(lWithoutVal);
                    end
                    parentObj.locData.loc.(['siteID_' 'masterAvg'])(idxNewLocs(lWithoutVal)) = dcal.site.ID;
                    parentObj.locData.loc.class(idxNewLocs(lWithoutVal)) = dcal.site.ID;
                    parentObj.locData.loc.(['rank_' 'masterAvg'])(idxNewLocs(lWithoutVal)) = dcal.site.indList;

                    % if one loc appears in more than one site, assign the spatial info
                    % based on the distance to origin.
                    if sum(~lWithoutVal)>0
                        idxWithVal = idxNewLocs(~lWithoutVal);
                        if isfield(newlocs,'znm')
                            distOriSq = parentObj.locData.loc.(['xnmaligned_' 'masterAvg'])(idxWithVal).^2+parentObj.locData.loc.(['ynmaligned_' 'masterAvg'])(idxWithVal).^2+parentObj.locData.loc.(['znmaligned_' 'masterAvg'])(idxWithVal).^2;
                            distNewSq = newlocs.xnm(~lWithoutVal).^2+newlocs.ynm(~lWithoutVal).^2+newlocs.znm(~lWithoutVal).^2;
                        else
                            distOriSq = parentObj.locData.loc.(['xnmaligned_' 'masterAvg'])(idxWithVal).^2+parentObj.locData.loc.(['ynmaligned_' 'masterAvg'])(idxWithVal).^2;
                            distNewSq = newlocs.xnm(~lWithoutVal).^2+newlocs.ynm(~lWithoutVal).^2;
                        end
                        idxSubWithVal = find(~lWithoutVal);
                        lReplace = distNewSq>distOriSq;
                        parentObj.locData.loc.(['xnmaligned_' 'masterAvg'])(idxWithVal(lReplace)) = newlocs.xnm(idxSubWithVal(lReplace));
                        parentObj.locData.loc.(['ynmaligned_' 'masterAvg'])(idxWithVal(lReplace)) = newlocs.ynm(idxSubWithVal(lReplace));
                        if isfield(newlocs,'znm')
                            parentObj.locData.loc.(['znmaligned_' 'masterAvg'])(idxWithVal(lReplace)) = newlocs.znm(idxSubWithVal(lReplace));
                        end
                        parentObj.locData.loc.(['siteID_' 'masterAvg'])(idxWithVal(lReplace)) = dcal.site.ID;
                        parentObj.locData.loc.class(idxWithVal(lReplace)) = dcal.site.ID;
                        parentObj.locData.loc.(['rank_' 'masterAvg'])(idxWithVal(lReplace)) = dcal.site.indList;
                    end
                end
            end
        end
        
        function [idxNewLocs, newlocs] = modifyMaster(obj, varargin)
            % to-do: LocMoFitGUI_2 here has to be replaced with
            % obj.currentFitter
            
            p = inputParser;
            fitter = obj.data(obj.currentData).fitter.LocMoFitGUI_2;
            roiSize = obj.parentObj.getPar('se_siteroi');
            p.addParameter('spatialOffset', 0);
            p.addParameter('binCloseAng', 0);
            p.addParameter('spatialTrimXY', zeros(1,2));
            p.addParameter('distFromOrigin', 1000);
            p.addParameter('binRadius', 300);
            p.addParameter('scalingFactor', 1);
            p.addParameter('firstBin', false);
            p.addParameter('siteID',0)
            p.parse(varargin{:});
            results = p.Results;
            
            parentObj = obj.parentObj;
            %% alignment           
%             if isempty(obj.idxCurrentSite)
%                 error('No ID of the current site specified.')
%                 return
%             end
            locData = parentObj.locData;
%             if length(obj.idxCurrentSite)==2
%                 locData.filter({'rank_masterAvg'},[],'minmax', obj.idxCurrentSite);
%                 % Here I use 'size' instead of 'position' because the filtering
%                 % is not based on the main xnm and ynm but rather xnmaligned_masterAvg','ynmaligned_masterAvg'
%                 [newlocs,indNewLocs] = locData.getloc({'xnmaligned_masterAvg','ynmaligned_masterAvg', 'znmaligned_masterAvg','layer'},'grouping', 'ungrouped', 'layer',find(locData.getPar('sr_layerson')),'removeFilter',{'filenumber'}); % per ROI info.
%             else
            [newlocs,indNewLocs] = locData.getloc({'xnmaligned_masterAvg','ynmaligned_masterAvg', 'znmaligned_masterAvg','layer','rank_masterAvg'},'grouping', 'ungrouped', 'layer',find(locData.getPar('sr_layerson')),'removeFilter',{'filenumber','rank_masterAvg'}); % per ROI info.
            lLocs = ismember(newlocs.rank_masterAvg,results.siteID);
            newlocs = subsetStruct(newlocs,lLocs);
            
            idxNewLocs = find(indNewLocs);
            idxNewLocs = idxNewLocs(lLocs);
            % convert
            % load allParsArg of the specific site

            fitter.externalInfo.binRadius = results.binRadius;
            fitter.externalInfo.scalingFactor = results.scalingFactor;
            fitter.externalInfo.binCloseAng = results.binCloseAng;
            fitter.convertNow
            
            newlocs.xnm = newlocs.xnmaligned_masterAvg-1000;
            newlocs.ynm = newlocs.ynmaligned_masterAvg;
            newlocs.znm = newlocs.znmaligned_masterAvg;
            
            % Transform the sites
            newlocs = fitter.locsRegister(newlocs, [], 1);
            
            % Trim locs
            indKept = abs(newlocs.xnm)<(roiSize/2 - results.spatialTrimXY(1));
            indKept = indKept&abs(newlocs.ynm)<(roiSize/2 - results.spatialTrimXY(2));
            
            newlocs = subsetStruct(newlocs,indKept);
            idxNewLocs = idxNewLocs(indKept);
            
            newlocs.xnm = newlocs.xnm+results.spatialOffset+results.distFromOrigin;
        end
        
        function saveToField(obj, filedName, idxLocs, locs)
            parentObj = obj.parentObj;
            lWithoutVal = parentObj.locData.loc.(['xnmaligned_' filedName])(idxLocs) == 0;
            parentObj.locData.loc.(['xnmaligned_' filedName])(idxLocs(lWithoutVal)) = locs.xnm(lWithoutVal);
            parentObj.locData.loc.(['ynmaligned_' filedName])(idxLocs(lWithoutVal)) = locs.ynm(lWithoutVal);
            if isfield(locs,'znm')
                parentObj.locData.loc.(['znmaligned_' filedName])(idxLocs(lWithoutVal)) = locs.znm(lWithoutVal);
            end
        end
        
        function oneFrame = mkFrame(obj, locs)
            p.sgauss=2;%smoothing of volume
            p.cutoff=4; % for isosurface, at max(V(:))/cutoff
            p.cutoffdc=[3 3]; % for isosurface, at max(V(:))/cutoff
            p.pxSize = 5; % pixel size
            
            mx = [0 200];
            my = [-200 200];
            mz = [-200 200];
            
            x = locs.xnm-1000;
            y = locs.ynm;
            z = locs.znm;
            x(x<0) = -x(x<0);
            [v, xb, yb,zb] = myhist3(x, y, z, p.pxSize, mx,my,mz);
            vs=smooth3(v,'gaussian',[1 1 1]*(round(2*p.sgauss/2)*2+1),p.sgauss);
            co=max(vs(:))/p.cutoff;
            f = obj.parentObj.getPar('hFrame');
            if isempty(f)||~isvalid(f)
                f = figure('Name','frame');
                obj.parentObj.setPar('hFrame',f);
            end
            clf(f)
            ax = axes(f);
            fc = [1,.75,.65];
%             fc = 'red';
            hiso = patch(ax, isosurface(vs,co),'EdgeColor','none','FaceColor',fc);
            isonormals(vs,hiso);
            hold on
            hcap=patch(ax, isocaps(vs,co),'FaceColor','interp','EdgeColor','none');
            axis equal
            xlim([0 range(my)/p.pxSize])
            ylim([0 range(mx)/p.pxSize])
            zlim([0 range(mz)/p.pxSize])
            view(-35,30)
            colormap hot
%             xlim([0,size(vs,1)])            
            lightangle(45,30);
            lighting gouraud
            colormap hot;
            drawnow
            oneFrame = getframe(f);
            oneFrame = oneFrame.cdata;
        end
        
        
    end
    methods (Static)
        function IDnewFormat = obj.setIDFormat(varargin)
            IDnewFormat = setIDFormat(varargin{:});
        end
    end
end

function IDnewFormat = setIDFormat(ID, format)
    switch format
        case 'meaning'
            IDnewFormat = regexprep(ID, '^LocMoFitGUI\_(\d)\.', '');
            IDnewFormat = regexprep(IDnewFormat, 'm\d{1,2}\.', '');
            IDnewFormat = regexprep(IDnewFormat, '\.l\d{1,2}$', '');
        case 'short'
            IDnewFormat = regexprep(ID, '^LocMoFitGUI\_(\d)\.', 'fit$1_');
            IDnewFormat = regexprep(IDnewFormat, '(m\d{1,2})\.', '$1_');
            IDnewFormat = regexprep(IDnewFormat, '\.(l\d{1,2})$', '_$1');
    end
end


function [f,imcomp]=plotcoords(x,y,z,R,thcov,p,x2,y2,z2)
% figure(88)
% scatter3(x,y,z+R,[],z+R)
% axis equal
% [t,p,r]=cart2sph(x,y,z);
% figure(89)
% polarplot(p,r,'.')
dz=R*cos(thcov);
dz=R;
xrange=[-p.minx,p.winsize];
yrange=[-p.winsize,p.winsize];
zrange=[-p.winsize,p.winsize]-p.winsize/2;

xr=xrange(1)+p.pixelsize:p.pixelsize:xrange(2)-p.pixelsize;
yr=yrange(1)+p.pixelsize:p.pixelsize:yrange(2)-p.pixelsize;
zr=zrange(1)+p.pixelsize:p.pixelsize:zrange(2)-p.pixelsize;
V=myhist3(x,y,-(z+dz),p.pixelsize,xrange,yrange,zrange);
%increase sampling by adding other half as well:
Vo=myhist3(x,y,-(z+dz),p.pixelsize,sort(-xrange),yrange,zrange);
V=V+Vo(end:-1:1,:,:);
Vs=smooth3(V,'gaussian',[1 1 1]*(round(2*p.sgauss/2)*2+1),p.sgauss);


V2=myhist3(x2,y2,-(z2+dz),p.pixelsize,xrange,yrange,zrange);
Vo2=myhist3(x2,y2,-(z2+dz),p.pixelsize,sort(-xrange),yrange,zrange);
V2=V2+Vo2(end:-1:1,:,:);
Vs2=smooth3(V2,'gaussian',[1 1 1]*(round(2*p.sgauss/2)*2+1),p.sgauss2);

[X,Y,Z]=meshgrid(yr,xr,zr);

co=max(Vs(:))/p.cutoff;
codc=max(Vs(:))/p.cutoffdc(1);

f1=figure(90);
clf
hiso=patch(isosurface(X,Y,Z,Vs,co),'EdgeColor','none','FaceColor',[1,.75,.65]);
isonormals(X,Y,Z,Vs,hiso);

hcap=patch(isocaps(X,Y,Z,Vs,co),'FaceColor','interp','EdgeColor','none');
axis equal
view(-35,30)

xlim(yrange)
ylim(xrange)
zlim(zrange)

lightangle(45,30);
lighting gouraud
colormap hot;


f2=figure(91);
clf
hiso=patch(isosurface(X,Y,Z,Vs,codc),'EdgeColor','none','FaceColor',[1,.8,.75],'FaceAlpha',1);
isonormals(X,Y,Z,Vs,hiso);
hold on


hcap=patch(isocaps(X,Y,Z,Vs,co),'FaceColor','interp','EdgeColor','none');
colormap hot
if 0
indin=x2>0;
scatter3(y2(indin),x2(indin),-z2(indin)-dz,2,'k')
else
codc2=max(Vs2(:))/p.cutoffdc(end);
% codc2=.5;
hiso2=patch(isosurface(X,Y,Z,Vs2,codc2),'EdgeColor','none','FaceColor',[0,1,1]);
isonormals(X,Y,Z,Vs2,hiso2);
hcap2=patch(isocaps(X,Y,Z,Vs2,codc2),'FaceColor',[0,.7,1],'EdgeColor','none');
end

axis equal
view(-35,30)

xlim(yrange)
ylim(xrange)
zlim(zrange)

lightangle(45,30);
lighting gouraud
whitebg('black')
f2.Children(1).CLim=f1.Children(1).CLim;

%2D analysis

img=getpolarimage(x,y,z,R,p)*p.polarfactor(1);
img2=getpolarimage(x2,y2,z2,R,p)*p.polarfactor(2);
s=size(img);
imcomp=zeros(s(2),s(1),3);
% imcomp2=imcomp;imcomp3=imcomp;
imcomp(:,:,1)=img';
imcomp(:,:,2)=img2';
% imcomp2(:,:,1)=img';
% imcomp3(:,:,2)=img2';
% implot=horzcat(imcomp,imcomp2,imcomp3);
f3=figure(45);
imagesc(imcomp);
axis equal

f=f2; %select which figure is saved
% indg=r~=0;
% % [z,x]=pol2cart(p(indg),r(indg));
% [z,x,y]=sph2cart(t(indg),p(indg),r(indg));
% zm=z+mean(r);
% figure(87);
% plot(x,-zm,'.')
% % xlim([-200 200])
% % ylim([-200 200])
% % posp.x=-zm;posp.y=x;
% ps=5;
% edgesx=-200:ps:200;
% edgesz=-100:ps:300;
% srim=histcounts2(x,zm,edgesx,edgesz);
% [Y,X]=meshgrid(rangez(1)+pixelsz/2:pixelsz:rangez(2),ranger(1)+pixelsx/2:pixelsx:ranger(2));

% srim=histrender(posp,[-200 200], [-200 200], 5, 5);
% figure(89)
% imagesc(srim')
% axis('equal')

% % subplot(2,3,1)
% 
% imn=srim./(X'+pixelsx/4);
% % imn=srim';
% imn2=[imn(:,end) imn(:,end:-1:2) imn];

% figure(88);polarplot(th,rch,'.')
% rlim([0,max(rmean)*1.1])
drawnow;
% waitforbuttonpress
% pause(.2)

end