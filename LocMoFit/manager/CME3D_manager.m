classdef CME3D_manager < SMLMModelFit_manager
    properties
    end
    methods
        function obj = CME3D_manager(varargin)
%             obj@SMLMModelFit_manager(1,2)
            obj@SMLMModelFit_manager(varargin{:});
        end
        function dynamicRec(obj)
            % Performing dynamic reconstructction based on the master
            % average (masterAvg)
            % this is now based on the masterAvg
            
            parentObj = obj.parentObj;
            p = obj.handles.linkedDynamicRec.getAllParameters;
            lGood = obj.filtering&obj.useSites;
            rank_site = find(lGood);
            numOfSite_filtered = sum(lGood);
            
            settings = obj.dynamicRec_settings(p);
            
            % reset the rules
            fitter = obj.data(obj.currentData).fitter.(settings.mainFitter);
            oldRules = fitter.converterRules;
            fitter.converterRules = [];
            fitter.converterUserDefined = [];
            
            for k = 1:length(settings.converter)
                fitter.converter([], settings.converter(k).rule, settings.converter(k).target)
            end
            
            % Bin the data based on their close angle
            radius = obj.getVariable_allSites('LocMoFitGUI_2.m1.radius');
            zOffset = obj.getVariable_allSites('LocMoFitGUI_2.m1.zOffset');
            theta = obj.getVariable_allSites('LocMoFitGUI_2.m1.closeAngle');
            closeAng = obj.getVariable_allSites('LocMoFitGUI_2.m1.realCloseAngle');
%             linkageError = obj.getVariable('LocMoFitGUI_2.m1.variation');
            binSize = floor(numOfSite_filtered/settings.binNumber);
            lastSiteInBins = binSize*settings.binNumber;
            rank_sitesInBin = rank_site(1:lastSiteInBins);
            %% Binning
%             obj.numOfSites binNumber
            rsRadius = reshape(radius(rank_sitesInBin), binSize,[]);
            rsTheta = reshape(theta(rank_sitesInBin), binSize,[]);
            rsZOffset = reshape(zOffset(rank_sitesInBin), binSize,[]);
            rsCloseAng = reshape(closeAng(rank_sitesInBin), binSize,[]);
%             rsLinkageError = reshape(linkageError(1:lastSite), binSize,[]);
            binZOffset = median(rsZOffset,1);
            binTheta = median(rsTheta,1);
            binRadius = median(rsRadius,1);
            binCloseAng = median(rsCloseAng,1);
%             binLinkageError = median(rsLinkageError,1);
%             binRadius = binRadius + binLinkageError/2;

            scalingFactor = binRadius./settings.masterAvgR;
            scalingFactor = reshape(scalingFactor, 1, [])';
            
            clipboard('copy', sprintf(['closeAngle\t' num2str(binTheta, '%.2f\t') '\nradius\t' num2str(binRadius, '%.2f\t')]));
            
            % reset the fields
            parentObj.locData.loc.(['xnmaligned_' 'masterAvgMod']) = zeros(size(parentObj.locData.loc.xnm));
            parentObj.locData.loc.(['ynmaligned_' 'masterAvgMod']) = zeros(size(parentObj.locData.loc.xnm));
            parentObj.locData.loc.(['znmaligned_' 'masterAvgMod']) = zeros(size(parentObj.locData.loc.xnm));
            
            % Skip the first bins, register the rest of bins with
            % re-scaling
            for k = settings.lastBin_noRescale+1:settings.binNumber
                siteID_oneBin = rank_sitesInBin((k-1)*binSize+1:k*binSize);
                currentBin = k;
                [idxNewLocs, locs] = obj.modifyMaster('spatialOffset', currentBin * settings.distBetweenBins+1000,...
                    'spatialTrimXY', settings.spatialTrimXY,...
                    'binRadius', binRadius(k),...
                    'scalingFactor', scalingFactor(k),...
                    'binCloseAng', binCloseAng(k),...
                    'siteID', siteID_oneBin);
                obj.saveToField('masterAvgMod', idxNewLocs, locs);
            end
            fitter.converterRules = oldRules;
            
            % Register the first bin without re-scaling
            et = 0;
            for currentBin = 1:settings.lastBin_noRescale
                for q = 1:binSize
                    real_siteID = (currentBin-1)*binSize + q;
                    obj.idxCurrentSite = rank_sitesInBin(real_siteID);
                    tic
                    obj.registerSites('firstBin', true,'spatialOffset', currentBin * settings.distBetweenBins+1000, 'spatialTrimXY', settings.spatialTrimXY);
                    et = et+toc;
                    if et >10
                        parentObj.status(['Register sites: Site ' num2str(q) ' of the ' num2str(binSize) ' sites in bin ' num2str(currentBin) ' done.']);
                        drawnow
                        et = 0;
                    end
                end
            end
            
        end
        
        
        function masterAvg(obj, varargin)
            p = inputParser;
            p.parse(varargin{:})
            p = p.Results;
            
            % Settings are defined in dynamicReconstruction_settings.
            settings = obj.masterAvg_settings;
            parentObj = obj.parentObj;
            
            mainFitter = settings.mainFitter;
            
            fitter = obj.data(obj.currentData).fitter.(mainFitter);
            
            % reset the rules
            oldRules = fitter.converterRules;
            fitter.converterRules = [];
            fitter.converterUserDefined = [];
            
            for k = 1:length(settings.converter)
                fitter.converter([], settings.converter(k).rule, settings.converter(k).target)
            end
            
            % reset the fields
            parentObj.locData.loc.(['xnmaligned_' 'masterAvg']) = zeros(size(parentObj.locData.loc.xnm));
            parentObj.locData.loc.(['ynmaligned_' 'masterAvg']) = zeros(size(parentObj.locData.loc.xnm));
            parentObj.locData.loc.(['znmaligned_' 'masterAvg']) = zeros(size(parentObj.locData.loc.xnm));
            parentObj.locData.loc.(['siteID_' 'masterAvg']) = zeros(size(parentObj.locData.loc.xnm));
            parentObj.locData.loc.(['rank_' 'masterAvg']) = zeros(size(parentObj.locData.loc.xnm));
            parentObj.locData.loc.class = zeros(size(parentObj.locData.loc.xnm));
            
            et = 0;
            for k = 1:sum(obj.useSites)
                obj.idxCurrentSite = k;
                tic
                obj.registerSites;
                et = et+toc;
                if et >10
                    parentObj.status(['Register sites: Site ' num2str(k) ' of the ' num2str(sum(obj.useSites)) ' sites done.']);
                    drawnow
                    et = 0;
                end
            end
            
            fitter.converterRules = oldRules;
        end
        
        function mkMovie(obj, varargin)
            % Last edits: 26.04.2022
            % Log:
            %   26.04.2022: fix the issue of black frames after the first
            %   frame.
            
            % If failed, roll back to dc91109e0bc306d9fce7f2d1aef6af451653576b
            inp = inputParser();
            inp.addParameter('saveTo','');
            inp.addParameter('numberOfSitesWithNegCur',0);
            inp.parse(varargin{:});
            results = inp.Results;
            mkMovSettings = obj.mkMovie_settings;
            % Performing dynamic reconstructction based on the master
            % average (masterAvg)
            % this is now based on the masterAvg
            dymSettings = obj.dynamicRec_settings;
            parentObj = obj.parentObj;
            
            lastSite = sum(obj.useSites);
            
            % reset the rules
            fitter = obj.data(obj.currentData).fitter.(dymSettings.mainFitter);
            oldRules = fitter.converterRules;
            fitter.converterRules = [];
            fitter.converterUserDefined = [];
            
            for k = 1:length(dymSettings.converter)
                fitter.converter([], dymSettings.converter(k).rule, dymSettings.converter(k).target)
            end
            
            % Bin the data based on their close angle
            curvature = obj.getVariable('LocMoFitGUI_2.m1.curvature');
            radius = obj.getVariable('LocMoFitGUI_2.m1.radius');
            closeAng = obj.getVariable('LocMoFitGUI_2.m1.realCloseAngle');
            numOfSites = length(radius);
%             lFlat = curvature<=0;
%             lFlat
            %% Binning
%             obj.numOfSites binNumber
            winSize = mkMovSettings.winSize;
            stepSize = mkMovSettings.stepSize;
            winBegin = 1:stepSize:numOfSites;
            winEnd = winSize:stepSize:numOfSites;
            winBegin = winBegin(1:length(winEnd));
            
            binRadius = zeros(size(winBegin));
            binCloseAng = zeros(size(winBegin));
            for k = 1:length(winBegin)
                binRadius(k) = mean(radius(winBegin(k):winEnd(k)));
                binCloseAng(k) = mean(closeAng(winBegin(k):winEnd(k)));
            end
%             binRadius = movmean(radius, winSize,'Endpoints', 'discard');
%             binCloseAng = movmean(closeAng, winSize, 'Endpoints', 'discard');
            
            scalingFactor = binRadius./dymSettings.masterAvgR;
            scalingFactor = reshape(scalingFactor, 1, [])';
            
            % !!! here the 2 should be 1 later
            
            %Jonas: introduce stepSize parameter (eg. in mkMovSettings)
            %make sure the 'layer' is saved. 
            par = figure;
            indSites = find(obj.usedSites);
            for k = 1:length(binRadius)
                indSites_bin = indSites(winBegin(k):winEnd(k));
                obj.idxCurrentSite = indSites_bin;
                if k ~= 1
                    [~, locs] = obj.modifyMaster('spatialOffset', 0, 'spatialTrimXY', [0 0], 'binRadius', binRadius(k), 'scalingFactor', scalingFactor(k), 'binCloseAng', binCloseAng(k), 'siteID',indSites_bin);
                else
                    % Register the first bin again without re-scaling
                    for l = indSites_bin
                        obj.idxCurrentSite = l;
                        if l == indSites_bin(1)
                            locs = obj.registerSites('firstBin', true,'spatialOffset', 0);
                        else
                            locsN = obj.registerSites('firstBin', true,'spatialOffset', 0);
                            locs = fuseParticle(locs, locsN);
                        end
                    end
                    locs.layer = ones(size(locs.xnm));
                end
                drawnow
                % calculate the time stamps
                if k==1
                    timeStamp = median([zeros(1, results.numberOfSitesWithNegCur) 1:(winSize-results.numberOfSitesWithNegCur)])/(numOfSites-results.numberOfSitesWithNegCur);
                    timeStamp = ['t = ' num2str(timeStamp,2)];
                    
                else
                    timeStamp = (mean([winBegin(k) winEnd(k)])-results.numberOfSitesWithNegCur)/(numOfSites-results.numberOfSitesWithNegCur);
                    timeStamp = ['t = ' num2str(timeStamp,2)];
                end
                oneFrame = render3DFrame(locs,par,'labelText',timeStamp);
                
                oneFrame = oneFrame.cdata;
                if k == 1
                    frames = uint8(zeros([size(oneFrame) length(binRadius)]));
                    frames(:,:,:,k) = oneFrame;
                else
                    frames(:,:,:,k) = oneFrame;
                end
                imwrite(frames(:,:,:,k), results.saveTo, 'Compression', 'none','WriteMode', "append");
                parentObj.status(['Make movies: Frame ' num2str(k) ' of the ' num2str(length(binRadius)) ' done.']);
            end
            implay(frames)         
            fitter.converterRules = oldRules;
            delete(par)
        end
        %%
        
        %% all settings
        function settings = masterAvg_settings(obj)
            % Here defines the converter for dynamic reconstruction
            % For generating master average
            
            settings.mainFitter = 'LocMoFitGUI_2';

            settings.converter(1).rule = 'find.m1.zOffset';
            settings.converter(1).target = 'post_z';
            
            settings.converter(end+1).rule = '150/find.m1.radius';
            settings.converter(end).target = 'post_scale';
        end
        
        function settings = firstBin_settings(obj)
            % Here defines the converter for dynamic reconstruction
            % For generating the first bin.
            
            settings.mainFitter = 'LocMoFitGUI_2';
            
            settings.converter(1).rule = 'find.m1.zOffset+find.m1.radius*sin(deg2rad(find.m1.realCloseAngle-90))';
            settings.converter(1).target = 'post_z';
        end
        
        function settings = dynamicRec_settings(obj,varargin)
            if isempty(varargin)
                % Here defines the converter for dynamic reconstruction
                % For dynamic reconstruction
                settings.binNumber = 10;
                settings.distBetweenBins = 350;
                settings.spatialTrimXY = [50 50];
                settings.masterAvgR = 150;

                settings.mainFitter = 'LocMoFitGUI_2';

    %             settings.converter(1).rule = '0';
    %             settings.converter(1).target = 'post_z';

                settings.converter(1).rule = 'find.scalingFactor';
                settings.converter(1).target = 'post_scale';
    %             
                settings.converter(end+1).rule = 'find.scalingFactor*(150*sin(deg2rad(find.binCloseAng)))';
    %             settings.converter(end+1).rule = 'find.scalingFactor*((-find.binRadius)*(1-sin(deg2rad(find.binCloseAng))))';            
                settings.converter(end).target = 'post_z';
            else
                p = varargin{1};
                % Here defines the converter for dynamic reconstruction
                % For dynamic reconstruction
                settings.binNumber = p.binNumber;
                settings.distBetweenBins = p.distBetweenBins;
                settings.spatialTrimXY = p.spatialTrimXY;
                settings.masterAvgR = p.masterAvgR;
                settings.lastBin_noRescale = p.lastBin_noRescale;
                settings.mainFitter = 'LocMoFitGUI_2';

    %             settings.converter(1).rule = '0';
    %             settings.converter(1).target = 'post_z';
                rules = p.convertTable.Data;
                for k = size(rules,1):-1:1
                    settings.converter(k).rule = rules{k,1};
                    settings.converter(k).target = rules{k,2};
                end
            end
        end
        
        function settings = mkMovie_settings(obj)
            % For making movie
            if isempty(obj.parentObj)
                settings.winSize = 30;
                settings.stepSize = 30;
            else
                p = obj.parentObj.getGuiParameters;
                settings.winSize = p.winSize;
                settings.stepSize = p.stepSize;
            end
            
        end
        
        
    end
end
