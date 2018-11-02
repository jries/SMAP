classdef CMErecWf < handle
    properties
        sourceType
        workFlow
        featureSet % might be called "originalFeatures" in the futuer
        featureMatrix
        featureRef
        dataSource
        temporalInfo % might be called "derivedFeatures" in the futuer
    end
    methods
        function obj = CMErecWf(se, isSimulation)
            obj.dataSource = se;
            if isSimulation
                obj.sourceType = 'simulation';

                numberOfSitesPerCell = obj.dataSource.numberOfSites/obj.dataSource.numberOfCells;
                siteIdx = 1:numberOfSitesPerCell:obj.dataSource.numberOfSites;
                [siteIdx,~] = meshgrid(siteIdx, 1:numberOfSitesPerCell);
                xData = siteIdx(:)';
                
                obj.setFeature('assignedTime', xData, 'Assigned time', 'Time (s)')
            else
                obj.sourceType = 'experiment';
            end
                        
            %% extract basic statistics
            
            sites=se().sites;
            numSites = se.numberOfSites;
            numSitesEachTP = se.numberOfSites/se.numberOfCells;
            siteIdx = 1:numSitesEachTP:numSites;
            [siteIdx,~] = meshgrid(siteIdx, 1:numSitesEachTP);
            siteIdx = siteIdx(:)';
            x1q13all = zeros(1, numSites);
            x2q13all = [];
            z1q13all = [];
            z2q13all = [];
            zmeddall = [];
            dPeakAll = zeros(1, numSites);
            N1all = zeros(1, numSites);
            N2all = zeros(1, numSites);
            disHat = [];
            for k=1:numSites
                try
                    x1q13all(k)=sites(k).evaluation.CME2CSide_yule2.x1q13;
                    x2q13=sites(k).evaluation.CME2CSide_yule2.x2q13;
                    z1q13=sites(k).evaluation.CME2CSide_yule2.z1q13;
                    z2q13=sites(k).evaluation.CME2CSide_yule2.z2q13;
                    zmedd=sites(k).evaluation.CME2CSide_yule2.zmd;
                    dPeak = max(sites(k).evaluation.CME2CSide_yule2.dPeak);
                    if length(dPeak)==1
                        dPeakAll(k) = dPeak;
                    end
                    N1all(k)=sites(k).evaluation.CME2CSide_yule2.N1;
                    N2all(k)=sites(k).evaluation.CME2CSide_yule2.N2;
                catch err
                    continue
                end
                x2q13all = [x2q13all; x2q13];
                z1q13all = [z1q13all; z1q13];
                z2q13all = [z2q13all; z2q13];
                zmeddall = [zmeddall; zmedd];
                
                % group locs of two proteins as one population, and get the
                % distance between to percentiles. This is an approximation of the range of proteins distributed in z
                
                % singleOneZ = [sites(k).evaluation.CME2CSide_yule2.locs1.ynmrot;sites(k).evaluation.CME2CSide_yule2.locs2.ynmrot];
                singleOneZ = [sites(k).evaluation.CME2CSide_yule2.z1;sites(k).evaluation.CME2CSide_yule2.z2];

                disHat(k) = prctile(singleOneZ, 99) - prctile(singleOneZ, 20);
            end
            %% set basic features
            obj.featureSet = [];
            obj.setFeature('x2q13all', x2q13all, 'Interquartile range (x) of protein2', 'nm')
            obj.setFeature('z1q13all', z1q13all, 'Interquartile range (z) of protein1', 'nm')
            obj.setFeature('z2q13all', z2q13all, 'Interquartile range (z) of protein2', 'nm')
            obj.setFeature('disHat', disHat', 'Distance (Prc90-Prc10)', 'nm')
            obj.setFeature('numberOfLocs', N1all'+N2all', 'Number of molecule', '')
            
        end
        
        function setFeature(obj, entity, data, label, unit)
            obj.featureSet.(entity) = [];
            obj.featureSet.(entity).data = data;
            obj.featureSet.(entity).label = label;
            obj.featureSet.(entity).unit = unit;
            
            [~,Idx] = sort(data);
            obj.featureSet.(entity).rank = Idx;
            
        end
        
        function rmFeature(obj, entity)
            obj.featureSet = rmfield(obj.featureSet, entity);
        end
        
        function setFeatureMatrix(obj, featureNames, matrixName)
            allFeatureNames=fieldnames(obj.featureSet);
            idx = find(ismember(allFeatureNames, featureNames));
            featureMatrixSingle = [];
            for i=1:size(idx)
                featureMatrixSingle = [featureMatrixSingle obj.featureSet.(allFeatureNames{idx(i)}).data];
            end
            obj.featureMatrix.(matrixName).source = featureNames;
            obj.featureMatrix.(matrixName).data = featureMatrixSingle;
        end
        
        function getRankIdx(obj, featureName) % TODO: will be removed soon
            [~,Idx] = sort(obj.featureSet.(featureName).data);
            obj.featureSet.(['rank_' featureName]) = [];
            obj.featureSet.(['rank_' featureName]).data = Idx;
            obj.featureSet.(['rank_' featureName]).unit = 'rank';
            obj.featureSet.(['rank_' featureName]).label = [obj.featureSet.(featureName).label ' (ranked)'];
        end
        
        function generateKdeImage(obj, imageSize, rangeOfSites, featureName, labelName, varargin)
            dataType = obj.whereFeature(featureName);
            idx = obj.(dataType).(featureName).rank;
            xnmrot1G = []; ynmrot1G = []; xnmrot2G = []; ynmrot2G = [];
            for i = rangeOfSites
                % for simulation
%                 xnmrot1G = [xnmrot1G; obj.dataSource.sites(idx(i)).evaluation.CME2CSide_yule2.locs1.xnmrot];
%                 ynmrot1G = [ynmrot1G; obj.dataSource.sites(idx(i)).evaluation.CME2CSide_yule2.locs1.ynmrot];
%                 xnmrot2G = [xnmrot2G; obj.dataSource.sites(idx(i)).evaluation.CME2CSide_yule2.locs2.xnmrot];
%                 ynmrot2G = [ynmrot2G; obj.dataSource.sites(idx(i)).evaluation.CME2CSide_yule2.locs2.ynmrot];
                
                xnmrot1G = [xnmrot1G; obj.dataSource.sites(idx(i)).evaluation.CME2CSide_yule2.x1];
                ynmrot1G = [ynmrot1G; obj.dataSource.sites(idx(i)).evaluation.CME2CSide_yule2.z1];
                xnmrot2G = [xnmrot2G; obj.dataSource.sites(idx(i)).evaluation.CME2CSide_yule2.x2];
                ynmrot2G = [ynmrot2G; obj.dataSource.sites(idx(i)).evaluation.CME2CSide_yule2.z2];
            end
            nm = length(xnmrot1G)+length(xnmrot2G);
            [groupDMap, bw]= getKernelMatrix(true, imageSize, [xnmrot1G ynmrot1G; xnmrot2G ynmrot2G], varargin{1,1}, varargin{1,2}); % Todo: the function to generate the image can be changed
            obj.featureRef.([featureName '_' labelName]) = [];
            obj.featureRef.([featureName '_' labelName]).nm = nm;
            obj.featureRef.([featureName '_' labelName]).Bandwidth = varargin{1,2};
            obj.featureRef.([featureName '_' labelName]).densityMap = groupDMap;
        end
        
        function measureSimilarity(obj, refName, imageSize, simFunc, coefBw, modeKdeMap) % TODO: customized imageSize
            %% Defined the way to define similarity using simFunc, and calculate similarity based on that.
            sim = [];
            for i = 1:obj.dataSource.numberOfSites
                Anm = obj.featureRef.(refName).nm;
                % for simulation
%                Bnm = length(obj.dataSource.sites(i).evaluation.CME2CSide_yule2.locs1.xnmrot) + length(obj.dataSource.sites(i).evaluation.CME2CSide_yule2.locs2.xnmrot);
                Bnm = length(obj.dataSource.sites(i).evaluation.CME2CSide_yule2.x1) + length(obj.dataSource.sites(i).evaluation.CME2CSide_yule2.x2);

                switch coefBw
                    case 0
                        cf = (Anm/Bnm)^0.5;
                        bw = cf*obj.featureRef.(refName).Bandwidth*2/5;
                    case 1
                        cf = (Anm/Bnm)^0.5;
                        bw = cf*obj.featureRef.(refName).Bandwidth*2;
                    case 2
                        cf = (Anm/Bnm)^0.5;
                        bw = cf*obj.featureRef.(refName).Bandwidth;
                    case 3
                        cf = (Anm/Bnm)^0.5;
                        bw = cf*obj.featureRef.(refName).Bandwidth*1/2;
                    case 4
                        cf = (Anm/Bnm)^0.5;
                        bw = cf*obj.featureRef.(refName).Bandwidth*1/3;
                    case 5
                        cf = (Anm/Bnm)^0.5;
                        bw = cf*obj.featureRef.(refName).Bandwidth*1/4;
                    case 6
                        cf = (Anm/Bnm)^0.5;
                        bw = cf*obj.featureRef.(refName).Bandwidth*1/5;
                    case 7
                        cf = (Anm/Bnm)^0.5;
                        bw = cf*obj.featureRef.(refName).Bandwidth*1/7;
                    case 8
                        cf = (Anm/Bnm)^0.5;
                        bw = cf*obj.featureRef.(refName).Bandwidth*1/8;
                    case 9
                        cf = (Anm/Bnm)^0.5;
                        bw = cf*obj.featureRef.(refName).Bandwidth*1/9;
                end
                % for simulation
%                locsAll = [obj.dataSource.sites(i).evaluation.CME2CSide_yule2.locs1.xnmrot obj.dataSource.sites(i).evaluation.CME2CSide_yule2.locs1.ynmrot; obj.dataSource.sites(i).evaluation.CME2CSide_yule2.locs2.xnmrot obj.dataSource.sites(i).evaluation.CME2CSide_yule2.locs2.ynmrot];
                locsAll = [obj.dataSource.sites(i).evaluation.CME2CSide_yule2.x1 obj.dataSource.sites(i).evaluation.CME2CSide_yule2.z1; obj.dataSource.sites(i).evaluation.CME2CSide_yule2.x2 obj.dataSource.sites(i).evaluation.CME2CSide_yule2.z2];
                switch modeKdeMap
                    case 'repeat'
                        
                        for ii = 1:20
                            subIdx = randsample(Bnm, round(Bnm*3/5));
                            locsSub = locsAll(subIdx,:);
                            pos = [];
                            pos.x = locsSub(:,1);
                            pos.y = locsSub(:,2);
                            pos.sx = zeros(size(pos.x),'like',pos.x);
                            pos.sy = zeros(size(pos.x),'like',pos.y);
                            pos.sx = (pos.sx + bw(1)).*1.5; 
                            pos.sy = (pos.sy + bw(2)).*1.5; 
                            pos.c = 0;
                            if ii == 1
                               singleOneSub = gaussrender_ellipt(pos, [-250 251], [-250 251], 1, 1 ,0,[0 0], 25);
                            else
                               singleOneSub(:,:,ii) = gaussrender_ellipt(pos, [-250 251], [-250 251], 1, 1 ,0,[0 0], 25);
                            end
                        end
                        singleOne = mean(singleOneSub, 3);
                        %figure(33354); imagesc(singleOne);
                        
                        
                    case 'repeat2'

                        for ii = 1:50
                            subIdx = randsample(Bnm, round(Bnm*1/2));
                            locsSub = locsAll(subIdx,:);
                            pos = [];
                            pos.x = locsSub(:,1);
                            pos.y = locsSub(:,2);
                            pos.sx = zeros(size(pos.x),'like',pos.x);
                            pos.sy = zeros(size(pos.x),'like',pos.y);
                            pos.sx = (pos.sx + bw(1)).*2; 
                            pos.sy = (pos.sy + bw(2)).*2; 
                            pos.c = 0;
                            if ii == 1
                               singleOneSub = gaussrender_ellipt(pos, [-250 251], [-250 251], 1, 1 ,0,[0 0], 25);
                            else
                               singleOneSub(:,:,ii) = gaussrender_ellipt(pos, [-250 251], [-250 251], 1, 1 ,0,[0 0], 25);
                            end
                        end
                        singleOne = mean(singleOneSub, 3);
                        figure(33354); imagesc(singleOne);
                    
                    case 'repeat3'
                        
                        for ii = 1:50
                            subIdx = randsample(Bnm, round(Bnm*4/5));
                            locsSub = locsAll(subIdx,:);
                            pos = [];
                            pos.x = locsSub(:,1);
                            pos.y = locsSub(:,2);
                            pos.sx = zeros(size(pos.x),'like',pos.x);
                            pos.sy = zeros(size(pos.x),'like',pos.y);
                            pos.sx = (pos.sx + bw(1)).*1.2; 
                            pos.sy = (pos.sy + bw(2)).*1.2; 
                            pos.c = 0;
                            if ii == 1
                               singleOneSub = gaussrender_ellipt(pos, [-250 251], [-250 251], 1, 1 ,0,[0 0], 25);
                            else
                               singleOneSub(:,:,ii) = gaussrender_ellipt(pos, [-250 251], [-250 251], 1, 1 ,0,[0 0], 25);
                            end
                        end
                        singleOne = mean(singleOneSub, 3);
                        figure(33354); imagesc(singleOne);
                        
                    case 'repeat4'
                        
                        for ii = 1:50
                            subIdx = randsample(Bnm, round(Bnm*3/5));
                            locsSub = locsAll(subIdx,:);
                            pos = [];
                            pos.x = locsSub(:,1);
                            pos.y = locsSub(:,2);
                            pos.sx = zeros(size(pos.x),'like',pos.x);
                            pos.sy = zeros(size(pos.x),'like',pos.y);
                            pos.sx = (pos.sx + bw(1)*1.1); 
                            pos.sy = (pos.sy + bw(2)*1.1); 
                            pos.c = 0;
                            if ii == 1
                               singleOneSub = gaussrender_ellipt(pos, [-250 251], [-250 251], 1, 1 ,0,[0 0], 25);
                            else
                               singleOneSub(:,:,ii) = gaussrender_ellipt(pos, [-250 251], [-250 251], 1, 1 ,0,[0 0], 25);
                            end
                        end
                        singleOne = mean(singleOneSub, 3);
                        figure(33354); imagesc(singleOne);
                    
                 	case 'normal'
                        pos = [];
                        pos.x = locsAll(:,1);
                        pos.y = locsAll(:,2);
                        pos.sx = zeros(size(pos.x),'like',pos.x);
                        pos.sy = zeros(size(pos.x),'like',pos.y);
                        pos.sx = pos.sx + bw(1); 
                        pos.sy = pos.sy + bw(2); 
                        pos.c = 0;
                    	singleOne = gaussrender_ellipt(pos, [-250 251], [-250 251], 1, 1 ,0,[0 0], 25);
                end

                %D = normxcorr2(avg, singleOne);
                %Ds = normxcorr2(SG, singleOne);
                %De = normxcorr2(EG, singleOne);
                %Dsa = normxcorr2(SGA, singleOne);
                %Dea = normxcorr2(EGA, singleOne);
                %Y(i) = max(D(:));
                %Ys(i) = max(Ds(:));
                %Ye(i) = max(De(:));
                %Ysa(i) = max(Dsa(:));
                %Yea(i) = max(Dea(:));

                switch simFunc
                case 'get3Dcorrshift'
                    [~,sim(i)]=get3Dcorrshift(obj.featureRef.(refName).densityMap, singleOne, 'max');
                case 'normxcorr2'
                    D=normxcorr2(obj.featureRef.(refName).densityMap, singleOne);
                    sim(i) = max(D(:));
                end
                componetsOfRefName = strsplit(refName,'_');
                obj.setFeature([refName '_Sim_bw' num2str(coefBw) '_mode_' modeKdeMap], sim', '', ['XCorr against ' componetsOfRefName{2}]);
            end
        end
        
        function temporalRecon(obj, featureMatrixName, methodName, basisOfS, zscored)
            featureMatrix = obj.featureMatrix.(featureMatrixName).data;
            dataType = obj.whereFeature(basisOfS);
            if zscored
                featureMatrix = zscore(featureMatrix,[],1);
            end
            switch methodName
                case 'wanderlust'
                    [~,sSite] = sort(obj.(dataType).(basisOfS).data);
                    sSite = sSite == 3; % the 3rd ranked site
                    wPar = [];
                    wPar.s = find(sSite);
                    rawWanderlust = wanderlust(featureMatrix, wPar);
                    temporalInfo = mean(rawWanderlust.traj);
            end
            if ~isempty(obj.temporalInfo)
                i = length(fieldnames(obj.temporalInfo))+1;
            else
                i = 1;
            end
            
            obj.temporalInfo.(['t' num2str(i)])=[]; % TODO: get number from metadata
            obj.temporalInfo.(['t' num2str(i)]).source = featureMatrixName;
            obj.temporalInfo.(['t' num2str(i)]).featureList = obj.featureMatrix.(featureMatrixName).source;
            obj.temporalInfo.(['t' num2str(i)]).method = methodName;
            obj.temporalInfo.(['t' num2str(i)]).data = temporalInfo;
            obj.temporalInfo.(['t' num2str(i)]).label = [featureMatrixName ' (' methodName ')'];
            obj.temporalInfo.(['t' num2str(i)]).unit = 'a.u.';
           
            [~,Idx] = sort(temporalInfo);
            obj.temporalInfo.(['t' num2str(i)]).rank = Idx;
        end
        
        function attr = whereFeature(obj, featureName)
            isFeature = isfield(obj.featureSet, featureName);
            isTempInfo = isfield(obj.temporalInfo, featureName);
            if isFeature
                attr = 'featureSet';
            elseif isTempInfo
                attr = 'temporalInfo';
            else
                attr = '';
            end
        end
            
        
        function getScatterPlot(obj, xName, yName)
            siteIdx = 1:obj.dataSource.numberOfSites;
            
            dataTypeX = obj.whereFeature(xName);
            dataTypeY = obj.whereFeature(yName);
           
            xTitle = [];
            if ~isequal(dataTypeX,'')
                xData = obj.(dataTypeX).(xName).data;
                xlabelName = obj.(dataTypeX).(xName).unit;
                xTitle = ['X: ' obj.(dataTypeX).(xName).label ' / Y: '];
                yData = obj.(dataTypeY).(yName).data;
            elseif isequal(xName,'rank')
                yData = obj.(dataTypeY).(yName).data;
                yData = yData(obj.(dataTypeY).(yName).rank);
                xlabelName = 'Rank';
                xData = siteIdx;
            end
            fullTitle = [xTitle obj.(dataTypeY).(yName).label];
            scatter(xData, yData);
            title(fullTitle)
            ylabel(obj.(dataTypeY).(yName).unit)
            xlabel(xlabelName)
        end

        function h = plotFeatureMatrix(obj, method, featureMatrixName, featureName, useRank)
            if useRank
               CS = obj.featureSet.(featureName).rank;
            else
               CS = obj.featureSet.(featureName);
            end
            
            switch method
                case 'tSNE'
                    h = tsne(obj.featureMatrix.(featureMatrixName));
                    figure; gscatter(h(:,1),h(:,2),CS, '','',10);

                case 'PCA'
                    [~,h,~] = pca(featureMatrixJCom);
                    figure; gscatter(h(:,1),h(:,2),CS, '','',10, 'off');
            end
        end
        
        function dataOutput = getData(obj, dataName, dataType)
            dataOutput = obj.(dataType).(dataName);
        end

    end
end


function [z, bw] = getKernelMatrix(centralized, matrixSize, sVrot, varargin) 
    [meshx, meshy] = meshgrid(0:matrixSize(1), 0:matrixSize(2)); % make full grid
    if centralized
        shx = matrixSize(1)/2; shy = matrixSize(2)/2;
    else
        shx = 0; shy = 0;
    end
    
    q = [meshx(:)-shx, meshy(:)-shy]; % convert the grid to positions
    
    if length(varargin)==0
        [sVrotx,xy,bw]=ksdensity(sVrot, q); % get kernel density estimation
    else
        [sVrotx,xy,bw]=ksdensity(sVrot, q, varargin{1,1}, varargin{1,2}); % get kernel density estimation    
    end
    
    % converted into an image
    idx = sub2ind(matrixSize(end:-1:1)+1, xy(:,2)+1+shx, xy(:,1)+1+shy); 
    z = zeros(matrixSize(end:-1:1)+1);
    z(idx) = sVrotx;
end