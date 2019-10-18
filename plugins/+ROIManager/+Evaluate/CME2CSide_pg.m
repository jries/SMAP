classdef CME2CSide_pg<interfaces.SEEvaluationProcessor
    properties
        savedevals
    end
    methods
        function obj=CME2CSide_pg(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function makeGui(obj)
            makeGui@interfaces.SEEvaluationProcessor(obj);
%             obj.guihandles.saveimagesb.Callback={@saveimagesb_callback,obj};
        end
        function out=run(obj,p)
            out=runintern(obj,p);
         
        end
        
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end

function pard=guidef
pard.reset.object = struct('Style','checkbox','String','Reset pos','Value',0);
pard.reset.position = [1,1];
pard.reset.Width = 1.5;

pard.doRealignment.object = struct('Style','checkbox','String','Re-alignment','Value',0);
pard.doRealignment.position = [2,1];
pard.doRealignment.Width = 1.5;

pard.t2.object = struct('Style','text','String','Z Ref point');
pard.t2.position = [2,2.5];
pard.t2.Width = 1;

pard.setMargin.object = struct('Style','edit','String','30');
pard.setMargin.position = [2,3.5];
pard.setMargin.Width = 0.5;

pard.useMeanShift.object = struct('Style','checkbox','String','Mean shift','Value',0);
pard.useMeanShift.position = [3,1];
pard.useMeanShift.Width = 1;

pard.doPlotMS.object = struct('Style','checkbox','String','Plot mean shift','Value',0);
pard.doPlotMS.position = [3,2];
pard.doPlotMS.Width = 1.5;

pard.plotCluster.object = struct('Style','checkbox','String','Plot clusters','Value',0);
pard.plotCluster.position = [3,3.5];
pard.plotCluster.Width = 1.5;

pard.doPlot.object = struct('Style','checkbox','String','Plot profile','Value',1);
pard.doPlot.position = [4,1];
pard.doPlot.Width = 1.5;

pard.t1.object = struct('Style','text','String','Profile cutoff (x,z)');
pard.t1.position = [5,1];
pard.t1.Width = 1.5;

pard.xrangeestDensityCutoff.object = struct('Style','edit','String','0.5');
pard.xrangeestDensityCutoff.position = [5,2.5];
pard.xrangeestDensityCutoff.Width = 0.5;

pard.yrangeestDensityCutoff.object = struct('Style','edit','String','0.75');
pard.yrangeestDensityCutoff.position = [5,3];
pard.yrangeestDensityCutoff.Width = 0.5;

pard.t3.object = struct('Style','text','String','Delt stop');
pard.t3.position = [6,1];
pard.t3.Width = 1;

pard.distStop.object = struct('Style','edit','String', '10');
pard.distStop.position = [6,2.5];
pard.distStop.Width = 0.5;

pard.t4.object = struct('Style','text','String','Max final shift');
pard.t4.position = [7,1];
pard.t4.Width = 1.5;

pard.finalShiftLimit.object = struct('Style','edit','String', '100');
pard.finalShiftLimit.position = [7,2.5];
pard.finalShiftLimit.Width = 0.5;

pard.plugininfo.type='ROI_Evaluate';

pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi','layer1_','layer2_','se_sitepixelsize'};
end



function out=runintern(obj,p)
%% Initiation

    p.doPlotCD=0;                       % show the cumulative density in Z 
    p.normalizeCD=0;                    % Normalize maximum to 1 
%     p.finalShiftLimit = 100;            % the maxima of the shift
    
%     p.doPlot = 1;                       % show the profile viewer
%     p.doPlotMS=0;                       % show the procedure of mean shift clustering
%     p.plotCluster = 0;                  % show the clusters identified
    
%     p.useMeanShift = 1;                 % enable mean shift
%     p.doRealignment = 0;                % enable realignment
    
%     p.setMargin=30;                     % set the margin
%     p.xrangeestDensityCutoff = 0.5;
%     p.yrangeestDensityCutoff = 0.75;
    
    out=[];
    roisize=ones(2,1)*p.se_siteroi;
    
    refyRoi = roisize(2)/2-p.setMargin; % Determine ref point of a ROI
    refxRoi = 0;
    modFlag = 0;                        
    
    if p.reset
        originalPos = obj.site.evaluation.CME2CSide_pg.originalPos;
        obj.site.pos = originalPos;
        out.originalPos = originalPos;
        return
    end
    
    if p.plotCluster                    %
        ff = figure(889878);
    end
    
    if ~isfield(obj.site.evaluation,'CME2CSide_pg')||~isfield(obj.site.evaluation.CME2CSide_pg,'originalPos')
        originalPos = obj.site.pos;
    else
        originalPos = obj.site.evaluation.CME2CSide_pg.originalPos;
    end
    
    while 1
        locs = obj.getLocs({'xnm','ynm','locprecnm','PSFxnm','xnmrot','ynmrot'},'group','grouped','layer',1,'size',obj.P.par.se_siteroi.content/2);
        
        if p.useMeanShift
            % find cluters
            [center, cluster,~] = MeanShiftCluster([locs.xnmrot locs.ynmrot]',50,p.doPlotMS);

            if p.plotCluster
                clf(ff, 'reset');
                axClu = axes(ff);
                scatter(axClu,locs.xnmrot,locs.ynmrot,[],cluster)
                axis(axClu, 'equal');
                set(axClu, 'Ydir', 'reverse');
            end

            ckept = eightCornerCheck(locs, center, cluster,15);
            usedMem = ismember(cluster, ckept);                     % locs in the cluster
            if p.plotCluster
                ffUsed = figure(9999);
                clf(ffUsed, 'reset');
                axUsed = axes(ffUsed);
                scatter(axUsed, locs.xnmrot(usedMem'), locs.ynmrot(usedMem'))
                axis(axUsed, 'equal')
                set(axUsed, 'Ydir', 'reverse')
            end
        else
            usedMem = logical(ones(size(locs.xnmrot))');
        end

        %histcounts()
        %profileViewer(locs, subset, nLocsLimit, binSize, setPad, p, doPlot, varargin)
        [siteProfile, qcFlag]= profileViewer(locs, usedMem', 25, 10, 30, p, p.doPlot, p.xrangeestDensityCutoff,p.yrangeestDensityCutoff);
        [siteProfileFull, ~]= profileViewer(locs, logical(ones(size(locs.xnmrot))), 25, 10, 30, p, 0, p.xrangeestDensityCutoff,p.yrangeestDensityCutoff);
        
        
        
        %% QC based on locs number
        if ~qcFlag
            obj.site.annotation.list2.value = 2;
            break
        else
            xRangeIndividual = findPeakRange(siteProfile.xdensity, siteProfile.xq, p.xrangeestDensityCutoff);
            xRangeEst = range(xRangeIndividual,2);
            if isempty(xRangeEst)
                xRangeEst = 0;
            end
            if any(xRangeEst> 200)
                obj.site.annotation.list2.value = 3;
                break
            end
            if max(xRangeEst)< 60
                obj.site.annotation.list2.value = 4;
                break
            end
        end
        
        %% find ref point
        % find x ref
%         factorD = 0.2;
%         maxyd = max(siteProfile.ydensity);
        cutoff = p.yrangeestDensityCutoff;
        peakRange = findPeakRange(siteProfile.ydensity, siteProfile.yq, cutoff);
        if isempty(peakRange)
            obj.site.annotation.list2.value = 7;
            break
        end
        midPointsY = movmean(siteProfile.yrange, 2);
        midPointsY = midPointsY(2:end);
        meaningfulRangeY = midPointsY>peakRange(1)&midPointsY<peakRange(2);
        meaningfulRangeY = meaningfulRangeY';
        meaningfulInd = zeros(size(meaningfulRangeY));
        meaningfulInd(siteProfile.grp) = meaningfulRangeY(siteProfile.grp);
        refxSite = mean(siteProfile.gm(ismember(siteProfile.grp, find(meaningfulInd))));
        
        % find y ref
%         refySite = prctile(locs.ynmrot(usedMem'), 95);
        refySite = max(peakRange(:));
        if isnan(refxSite)||isnan(refySite)
            obj.site.annotation.list2.value = 7;
            break
        end
        sxRoiRot = refxSite - refxRoi;
        syRoiRot = refySite - refyRoi;
        
        if ~p.doRealignment
            break
        end
        

        
        if norm([sxRoiRot syRoiRot; 0 0])<p.distStop
            break
        elseif modFlag>=20
            obj.site.annotation.list2.value = 5;
            break
        end
        
        %% update the center of the site
        obj.site.annotation.rotationpos.angle = pos2angle(obj.site.annotation.rotationpos.pos);
        [sxRoi, syRoi] = rotcoord(sxRoiRot, syRoiRot, (-(obj.site.annotation.rotationpos.angle)/180)*pi);

        newPos = obj.site.pos + [sxRoi syRoi 0];
        
        % if the shift of ROI center is too large, disable the ROI
        if norm(originalPos - newPos) > p.finalShiftLimit;
            obj.site.annotation.list2.value = 6;
            break
        end
        obj.site.pos = newPos;
        pause(10e-1000);
        modFlag = modFlag+1;
    end
    
    if qcFlag
    % estimate the radius of sites
        factorD = 0.2;
        maxxd = max(siteProfile.xdensity);
        cutoff = maxxd*factorD;
        peakRange = findPeakRange(siteProfile.xdensity, siteProfile.xq, cutoff);
        
        
        estDia = range(peakRange(:)); % measure the diameter (the name should be chaged: estRad->estDia) by the x profile
        
        estYR = range(prctile(locs.ynmrot(usedMem'), [15 95])); % estimate the z range based on major cluster(s)
%         aucDensity = sum(cumDensity(siteProfile, p));
        ypeaks = findPeakRange(siteProfile.ydensity, siteProfile.yq, p.yrangeestDensityCutoff);
        estYRbyDensity = range(ypeaks(:));
        
        % model fitting
%         ftx = fittype( 'xGauss(x, a1, b1, c1, d, aF, c3)', 'coefficients',{'a1', 'b1', 'c1', 'd', 'aF', 'c3'},'independent','x', 'dependent', 'y');
%         fx = fit( siteProfile.xq', siteProfile.xdensity', ftx, 'StartPoint', [0.01, -70, 35, 70, 0.02, 35] , 'Lower', [0, -150, 20,5, 1,2.5], 'Upper', [inf, 0, 150, 150, 10,inf]);
        %figure(9895); plot( fx, siteProfile.xq', siteProfile.xdensity' )
        
%         fty = fittype( 'yGauss(x, a1, b1, c1, d, aF, c2)', 'coefficients',{'a1', 'b1', 'c1', 'd', 'aF', 'c2'},'independent','x', 'dependent', 'y');
%         fy = fit( siteProfile.yq', siteProfile.ydensity', fty, 'StartPoint', [0.01, 100, 20, -70, 0.01, 20], 'Lower', [0, 0, 2.5,-300, 0,2.5], 'Upper', [inf, 150, 70, 0, 10,70]);
        %figure(99); plot( fy, siteProfile.yq', siteProfile.ydensity' )
        
        out.estimatedDiameter = estDia;
        out.estimatedZRange = estYR;
%         out.aucDensity = aucDensity;
        out.estimatedZRangeByDensity = estYRbyDensity;
%         out.fy = fy;
%         out.fx = fx;
        out.siteProfile = siteProfile;
        out.filter = usedMem;
        out.xRangeEst = xRangeEst;
    end
    out.originalPos = originalPos;
    %obj.site.image = [];
    %obj.locData.SE.plotsite(obj.site,obj.locData.SE.processors.preview.guihandles.siteax,[]);
    %drawnow
    %obj.locData.SE.processors.preview.updateSitelist;
end

function peakRange = findPeakRange(densityFunc, query, cutoff)
    aboveCut = densityFunc > cutoff;
    shiDiff = [aboveCut 0] - [0 aboveCut];
    indRise = find(shiDiff == 1);
    indFall = find(shiDiff == -1)-1;
    peakRange = [indRise' indFall'];
    peakRange = [query(peakRange(:,1))' query(peakRange(:,2))'];
end

function theRange = getRange(locs, what2get, binSize)
    if ~isempty(what2get)
        theRange = floor(min(locs.(what2get))/binSize)*binSize:binSize:ceil(max(locs.(what2get))/binSize)*binSize;
    else
        theRange = floor(min(locs)/binSize)*binSize:binSize:ceil(max(locs)/binSize)*binSize;
    end
end

function ckept = eightCornerCheck(locs, center, cluster, numLimit)
numOfCluster = size(center, 2);
numelCluster = grpstats(cluster,cluster,'numel');
numOfClusterBylabel = length(numelCluster);
NOCAgree= numOfClusterBylabel==numOfCluster ;
ckeptP1 = [];
ckept = 1;
    if numOfCluster > 1 && NOCAgree
        while 1
            numelCluster = grpstats(cluster,cluster,'numel');
            if (size(center, 2))==1
                break
            end
            [~,maxInd] = max(numelCluster);
            remainc = unique(cluster);
            maxInd = remainc(maxInd);
            % connected by x
            connected1 = rangesearch(center(1,maxInd)',center(1,[1:maxInd-1 maxInd+1:end])',70, 'NSMethod', 'exhaustive');
            % connected by z
            connected2 = rangesearch(center(2,maxInd)',center(2,[1:maxInd-1 maxInd+1:end])',20, 'NSMethod', 'exhaustive');
            connected = zeros([1 numOfCluster]);
            if (isempty(connected1{1})||isempty(connected2{1}))
                break
            else
                selected = [connected1{1} connected2{1}];
                selected(selected>=maxInd) = selected(selected>=maxInd)+1;
                for k = 1:length(selected)
                    connected(selected(k)) = connected(selected(k))+1;
                end
                connected = find(connected==2);
            end

            mergeFlag = [maxInd connected];
            newCenter = mean(numelCluster(mergeFlag) .* center(:,mergeFlag)')';

            cluster(cluster==connected)=maxInd;

            center(:,maxInd)=newCenter;
            center(:,connected)=[-300;-300];
            ckeptP1 = [ckeptP1 connected];
        end
        cxLower = center(1,maxInd)-70;
        cxUpper = center(1,maxInd)+70;
        ckept = find(center(1,:)>=cxLower&center(1,:)<=cxUpper);
    end
    
    if numLimit>0
        c2rm = find(numelCluster < numLimit);
        remainc = unique(cluster);
        c2rm = remainc(c2rm)';
        c2rmInd = ismember(ckept, c2rm);
        ckept = ckept(~c2rmInd);
    end
    ckept = unique([ckept ckeptP1]);
end

function [out,qcFlag] = profileViewer(locs, subset, nLocsLimit, binSize, setPad, p, doPlot, varargin)
    % yrange = getRange(locs, 'ynmrot', binSize);
    % xrange = getRange(locs, 'xnmrot', binSize);
    if sum(subset)<nLocsLimit
        qcFlag = 0;
        out = [];
    else
        qcFlag = 1;
        roisize=ones(2,1)*p.se_siteroi;
    
        xlocsKept = locs.xnmrot(subset');
        ylocsKept = locs.ynmrot(subset');
        nLocs = length(xlocsKept);
        xrange = getRange(xlocsKept, [], binSize);
        yrange = getRange(ylocsKept, [], binSize);

        [~,~,bin] = histcounts(ylocsKept ,yrange);

        [gm,grp]=grpstats(xlocsKept, bin, {'mean', 'gname'});
        grp = cellfun(@str2num, grp);
        [xdensity,xq]=ksdensity(xlocsKept,round(min(xlocsKept)-setPad):round(max(xlocsKept)+setPad), 'Bandwidth',10);
        [ydensity,yq]=ksdensity(ylocsKept,round(min(ylocsKept)-setPad):round(max(ylocsKept)+setPad), 'Bandwidth',10);

        out.xdensity = xdensity*nLocs; out.xq = xq; out.ydensity=ydensity*nLocs; out.yq=yq; out.gm=gm; out.grp=grp; out.bin=bin; out.xrange=xrange; out.yrange=yrange; out.padding = setPad;
        if doPlot
            viewPxSize = 2;
            % img = histcounts2(locs.ynmrot, locs.xnmrot, getRange(locs, 'ynmrot', viewPxSize), getRange(locs, 'ynmrot', viewPxSize));
            % h = fspecial('gaussian', 10, 10);
            % imgf = imfilter(img, h);
            f998 = figure(998); 
            % imagesc(imgf);
            % colormap(f998,mymakelut('red hot'))

            f998;
            % ax998_4=subplot(2,2,4);
            ax998_4=subplot('position',[0.3 0.05 0.65 0.65]);
            scatter(ax998_4,xlocsKept, ylocsKept)

    %         grp = cellfun(@str2num,grp);
            hold on
            shift = (roisize(1)/2)/viewPxSize;
            plot(ax998_4, [gm'; gm'],[yrange(grp); yrange(grp+1)])
            set(ax998_4, 'Ydir', 'reverse')
            axis(ax998_4, 'equal')
            hold off


            % ax998_2=subplot(2,2,2);
            ax998_2=subplot('position',[0.3 0.75 0.65 0.2]);
            
            if ~isempty(varargin)
                plot(ax998_2,xq,out.xdensity, '-b',xq,repelem(varargin{1}, length(xq)), '-r')
            else
                plot(ax998_2,xq,out.xdensity)
            end
            axis(ax998_2, [ax998_4.XLim ax998_2.YLim])
            

            % ax998_3=subplot(2,2,3);
            ax998_3=subplot('position',[0.05 0.05 0.2 0.65]);
            if ~isempty(varargin)
                plot(ax998_3,yq,out.ydensity, '-b',yq,repelem(varargin{2}, length(yq)), '-r')
            else
                plot(ax998_3,yq,out.ydensity)
            end
            set(ax998_3, 'Ydir', 'reverse')
            axis(ax998_3, [ax998_4.YLim ax998_3.YLim])
            
            view(ax998_3,90,90)
        end
    end
    
end

function CD = cumDensity(siteProfile,p)
    roiRange = p.se_siteroi+siteProfile.padding*6;
    CD = zeros([roiRange 1]);
    CD(-siteProfile.yq+roiRange/2) = siteProfile.ydensity;
    CD = cumsum(CD);
    if p.normalizeCD
        CD = CD/max(CD);
    end
    if p.doPlotCD
        figure(99989); 
        hold on
        plot((1:roiRange)', CD)
        hold off
    end
end
