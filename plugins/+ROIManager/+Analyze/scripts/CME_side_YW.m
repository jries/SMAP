global se
sites=se.sites;
if 0
    %% extract basic statistics
    numSites = length(sites);
    siteIdx = 1:30:numSites;
    [siteIdx,~] = meshgrid(siteIdx, 1:30);
    siteIdx = siteIdx(:)';
    x1q13all = zeros(1, numSites);
    x2q13all = [];
    z1q13all = [];
    z2q13all = [];
    zmeddall = [];
    dPeakAll = zeros(1, numSites);
    N1all = zeros(1, numSites);
    N2all = zeros(1, numSites);
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
    end

    %% the range of proteins distributed in z
    disHat = [];
    for i = 1:numSites
        singleOneZ = [sites(i).evaluation.CME2CSide_yule2.locs1.ynmrot;sites(i).evaluation.CME2CSide_yule2.locs2.ynmrot];
        disHat(i) = prctile(singleOneZ, 90) - prctile(singleOneZ, 10);
    end
    [~,IdisHat] = sort(disHat');

    %% Prepare references for comparison based on xcorr
    SG = calGruopRefDMap(IdisHat, 1:100);
    EG = calGruopRefDMap(IdisHat, numSites:-1:numSites-100);
    SGA = calGruopRefDMap(1:numSites, 1:100);
    EGA =  calGruopRefDMap(1:numSites, numSites:-1:numSites-100);
    MG = calGruopRefDMap(IdisHat, round(numSites/2)-49:round(numSites/2)+50);

    while 0
        figure(20400)
        opts = optimoptions('ga','PlotFcn',@gaplotbestf);
        [x, yval, fg]= ga(@findMinCor, 300, [repelem(1,300);repelem(-1,300)],[150;-150],[],[],repelem(0,300),repelem(1,300),[],1:300,opts);
        xy1G1 = getKernelMatrix(true, [500 500], [xnmrot1G1 ynmrot1G1]);
        xy2G1 = getKernelMatrix(true, [500 500], [xnmrot2G1 ynmrot2G1]);
        xy1G2 = getKernelMatrix(true, [500 500], [xnmrot1G2 ynmrot1G2]);
        xy2G2 = getKernelMatrix(true, [500 500], [xnmrot2G2 ynmrot2G2]);
        figure(2002)
        subplot(2,2,1);
        imagesc(xy1G1)
        subplot(2,2,2);
        imagesc(xy2G1)
        subplot(2,2,3);
        imagesc(xy1G2)
        subplot(2,2,4);
        imagesc(xy2G2)

        xy1D = xcorr2(xy1G1,xy1G2);
        max(xy1D(:))

        xy2D = xcorr2(xy2G1,xy2G2);
        max(xy2D(:))

    end

    %Y = [];
    %Ys = [];
    %Ye = [];
    %Ysa = [];
    %Yea = [];

    %% get the similarity based on xcross agains different references
    Js = []; Je = []; Jm = []; Jsa = []; Jea = [];
    for i = 1:numSites
        singleOne = getKernelMatrix(true, [500 500], [sites(i).evaluation.CME2CSide_yule2.locs1.xnmrot sites(i).evaluation.CME2CSide_yule2.locs1.ynmrot; sites(i).evaluation.CME2CSide_yule2.locs2.xnmrot sites(i).evaluation.CME2CSide_yule2.locs2.ynmrot]);
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
        [~,Js(i)]=get3Dcorrshift(SG, singleOne, 'max');
        [~,Je(i)]=get3Dcorrshift(EG, singleOne, 'max');
        [~,Jm(i)]=get3Dcorrshift(MG, singleOne, 'max');
        [~,Jsa(i)]=get3Dcorrshift(SGA, singleOne, 'max');
        [~,Jea(i)]=get3Dcorrshift(EGA, singleOne, 'max');
    end
end
CMErec = CMErecWf(se, true);
figure; CMErec.getScatterPlot('rank', 'disHat');
CMErec.generateKdeImage([500 500], 1:15, 'disHat', 'top15', 'Bandwidth', [10 10]);
CMErec.generateKdeImage([500 500], CMErec.dataSource.numberOfSites:-1:(CMErec.dataSource.numberOfSites-15), 'disHat', 'bottom15', 'Bandwidth', [10 10]);
CMErec.generateKdeImage([500 500], (round(CMErec.dataSource.numberOfSites/2)-7):(round(CMErec.dataSource.numberOfSites/2)+8), 'disHat', 'median15', 'Bandwidth', [10 10]);

figure; imagesc(CMErec.featureRef.disHat_top15.densityMap)
figure; imagesc(CMErec.featureRef.disHat_median15.densityMap)
figure; imagesc(CMErec.featureRef.disHat_bottom15.densityMap)


CMErec.measureSimilarity('disHat_top15', [500 500], 'normxcorr2', 3, 'normal')
CMErec.measureSimilarity('disHat_bottom15', [500 500], 'normxcorr2', 3, 'normal')
CMErec.measureSimilarity('disHat_median15', [500 500], 'normxcorr2', 3, 'normal')



figure; CMErec.getScatterPlot('disHat_top15_Sim_bw3_mode_normal', 'disHat_bottom15_Sim_bw3_mode_normal')
%CMErec.setFeatureMatrix(["disHat" "disHat_bottom30_Sim_bw3_mode_normal" "disHat_median30_Sim_bw3_mode_normal" "disHat_top30_Sim_bw3_mode_normal"], 'A')
CMErec.setFeatureMatrix(["disHat" "numberOfLocs" "disHat_bottom15_Sim_bw3_mode_normal" "disHat_median15_Sim_bw3_mode_normal" "disHat_top15_Sim_bw3_mode_normal"], 'B')

%CMErec.temporalRecon('A', 'wanderlust', true)
CMErec.temporalRecon('B', 'wanderlust', 'disHat', true)

figure; CMErec.getScatterPlot('rank','t1')
%figure; CMErec.getScatterPlot('assignedTime','t1')
figure; CMErec.getScatterPlot('disHat', 't1')
%figure; CMErec.getScatterPlot('assignedTime','disHat')

CMErec.generateKdeImage([500 500], 1:15, 't1', 'top15', 'Bandwidth', [10 10]);
CMErec.generateKdeImage([500 500], CMErec.dataSource.numberOfSites:-1:(CMErec.dataSource.numberOfSites-15), 't1', 'bottom15', 'Bandwidth', [10 10]);
CMErec.generateKdeImage([500 500], (round(CMErec.dataSource.numberOfSites/2)-7):(round(CMErec.dataSource.numberOfSites/2)+8), 't1', 'median15', 'Bandwidth', [10 10]);

figure; imagesc(CMErec.featureRef.t1_top15.densityMap)
figure; imagesc(CMErec.featureRef.t1_median15.densityMap)
figure; imagesc(CMErec.featureRef.t1_bottom15.densityMap)



for i = 1:se.numberOfSites
    se.sites(i).evaluation.CME2CSide_yule2.t1 = CMErec.temporalInfo.t1.data(i);
    se.sites(i).evaluation.CME2CSide_yule2.disHat = CMErec.featureSet.disHat.data(i);
end

% CMErec = CMErecWf(se);
% CMErec.getRankIdx('disHat');
% CMErec.getFeaturePlot('disHat', 'rank');
% CMErec.generateKdeImage([500 500], 1:30, 'disHat', 'top30', 'Bandwidth', [10 10]);
% CMErec.generateKdeImage([500 500], CMErec.dataSource.numberOfSites:-1:(CMErec.dataSource.numberOfSites-30), 'disHat', 'bottom30', 'Bandwidth', [10 10]);
% CMErec.generateKdeImage([500 500], (round(CMErec.dataSource.numberOfSites/2)-14):(round(CMErec.dataSource.numberOfSites/2)+15), 'disHat', 'median30', 'Bandwidth', [10 10]);
% for i=1:9
%     CMErec.measureSimilarity('disHat_top30', [500 500], 'normxcorr2', i, 'normal')
% end
% CMErec.measureSimilarity('disHat_median30', [500 500], 'normxcorr2', 3, 'normal')
% CMErec.measureSimilarity('disHat_bottom30', [500 500], 'normxcorr2', 3, 'normal')
% CMErec.measureSimilarity('disHat_top30', [500 500], 'normxcorr2', 4, 'repeat')
% CMErec.measureSimilarity('disHat_top30', [500 500], 'normxcorr2', 4, 'repeat2')
% figure; CMErec.getFeaturePlot('disHat_top30_Sim_bw4_mode_repeat', 'assignedTime')
% figure; CMErec.getFeaturePlot('disHat_top30_Sim_bw4_mode_normal', 'assignedTime')
% figure; CMErec.getFeaturePlot('disHat', 'assignedTime')
% CMErec.getRankIdx('disHat_top30_Sim_bw4_mode_repeat');
% CMErec.getRankIdx('disHat_top30_Sim_bw4_mode_normal');
% CMErec.getRankIdx('disHat_top30_Sim_bw4_mode_repeat2');
% figure; CMErec.getFeaturePlot('rank_disHat_top30_Sim_bw4_mode_repeat2', 'assignedTime')
% hold on
% CMErec.getScatterPlot('rank_disHat_top30_Sim_bw4_mode_repeat', 'assignedTime', 'featureSet')
% %CMErec.getFeaturePlot('rank_disHat_top30_Sim_bw4_mode_normal', 'assignedTime')
% hold off
% CMErec.setFeatureMatrix(["disHat" "disHat_bottom30_Sim_bw3_mode_normal" "disHat_median30_Sim_bw3_mode_normal" "disHat_top30_Sim_bw3_mode_normal"], 'A')
% CMErec.temporalRecon('A', 'wanderlust', true)
% CMErec.getScatterPlot('t1', 'assignedTime', 'temporalInfo')

if 0
    %% Create reference for the plots
    idxTime1 = se.P.par.proteinES.content.sla2.timePar.time<100;
    time1 = se.P.par.proteinES.content.sla2.timePar.time(idxTime1);
    nm1 = se.P.par.proteinES.content.sla2.timePar.nm(idxTime1);
    zm1 = se.P.par.proteinES.content.sla2.timePar.zm(idxTime1);
    de1 = se.P.par.proteinES.content.sla2.timePar.de(idxTime1);

    idxTime2 = se.P.par.proteinES.content.las17.timePar.time<100;
    time2 = se.P.par.proteinES.content.las17.timePar.time(idxTime2);
    nm2 = se.P.par.proteinES.content.las17.timePar.nm(idxTime2);
    zm2 = se.P.par.proteinES.content.las17.timePar.zm(idxTime2);
    de2 = se.P.par.proteinES.content.las17.timePar.de(idxTime2);

    refZMedDal = zm1-zm2;

    time1 = (time1)*30/0.2;
    time2 = (time2)*30/0.2;

    %% Try to check some features

    figure(500)
    subplot(2,5,1);
    scatter(siteIdx, N1all'+N2all');
    title('Number of localizations')
    xticks(time1)
    xticklabels(time1/30*0.2)
    hold on
    plot(time1,nm1+nm2,'LineWidth',2, 'Color', 'r')
    hold off

    subplot(2,5,2);
    scatter(siteIdx, disHat);
    title('Distance (Prc90-Prc10)')
    ylabel('nm')
    xticks(time1)
    xticklabels(time1/30*0.2)
    hold on
    plot(time2,refZMedDal,'LineWidth',2, 'Color', 'r')
    hold off

    subplot(2,5,3);
    scatter(siteIdx, z1q13all);
    ylabel('nm')
    title('Interquartile range (z) of protein1')
    xticks(time1)
    xticklabels(time1/30*0.2)
    hold on
    plot(time1,de1,'LineWidth',2, 'Color', 'r')
    hold off

    subplot(2,5,6);
    scatter(siteIdx, z2q13all);
    ylabel('nm')
    title('Interquartile range (z) of protein2')
    xticks(time1)
    xticklabels(time1/30*0.2)
    hold on
    plot(time2,de2,'LineWidth',2, 'Color', 'r')
    hold off

    subplot(2,5,7);
    scatter(siteIdx, Jm);
    title('XCorr against average')
    xticks(time1)
    xticklabels(time1/30*0.2)

    subplot(2,5,4);
    scatter(siteIdx, Js);
    title('JXCorr against top 100')
    xticks(time1)
    xticklabels(time1/30*0.2)

    subplot(2,5,5);
    scatter(siteIdx, Je);
    title('XCorr against bottom 100')
    xticks(time1)
    xticklabels(time1/30*0.2)

    subplot(2,5,9);
    scatter(siteIdx, Jsa);
    title('XCorr against topA 100')
    xticks(time1)
    xticklabels(time1/30*0.2)

    subplot(2,5,10);
    scatter(siteIdx, Jea);
    title('XCorr against bottomA 100')
    xticks(time1)
    xticklabels(time1/30*0.2)

    %% ranking based on different features
    %[~,IY] = sort(Y');
    %[~,IYs] = sort(Ys');
    %[~,IYe] = sort(Ye');
    NmAll = N1all'+N2all';
    %[~,INm] = sort(NmAll);

    %% generate feature matrics
    %featureMatrixCom = [Y' Ys' Ye' disHat' NmAll];
    featureMatrixJCom = [zscore(Jm') zscore(Js') zscore(Je') zscore(disHat') zscore(NmAll)];
    featureMatrixComI = [IY IYs IYe IdisHat INm];
    featureMatrixComZ = [zscore(Y)' zscore(Ys)' zscore(Ye)' zscore(disHat)' zscore(NmAll)];
    featureMatrixComZA = [zscore(Y)' zscore(Ysa)' zscore(Yea)' zscore(disHat)' zscore(NmAll)];
    featureMatrixXcorrZ = [zscore(Y)' zscore(Ys)' zscore(Ye)'];
    featureMatrixAbi1 = [zscore(Ys)' zscore(Ye)' zscore(disHat)' zscore(NmAll)];
    featureMatrixAbi2 = [zscore(Ys)' zscore(Ye)' zscore(disHat)'];

    %featureMatrix = [Y' Ys' Ye' disHat' N1all'+N2all'];

    figure(501);
    subplot(2,3,1)
    tSNECom = tsne(featureMatrixCom);
    gscatter(tSNECom(:,1),tSNECom(:,2),siteIdx, '','',10, 'off');
    title('SimA, SimS, SimE, disHat, & NmAll (Raw)')

    subplot(2,3,2)
    tSNEComI = tsne(featureMatrixComI);
    gscatter(tSNEComI(:,1),tSNEComI(:,2),siteIdx, '','',10, 'off');
    title('SimA, SimS, SimE, disHat, & NmAll (Ranked)')

    subplot(2,3,3)
    tSNEComZ = tsne(featureMatrixComZ);
    gscatter(tSNEComZ(:,1),tSNEComZ(:,2),siteIdx, '','',10, 'off');
    title('SimA, SimS, SimE, disHat, & NmAll (Z-scored)')

    subplot(2,3,4)
    tSNEXcorrZ = tsne(featureMatrixXcorrZ);
    gscatter(tSNEXcorrZ(:,1),tSNEXcorrZ(:,2),siteIdx, '','',10, 'off');
    title('SimA, SimS, & SimE (Z-scored)')

    subplot(2,3,5)
    tSNEAbi1= tsne(featureMatrixAbi1);
    gscatter(tSNEAbi1(:,1),tSNEAbi1(:,2),siteIdx, '','',10, 'off');
    title('SimS, SimE, disHat, & NmAll (Z-scored)')

    subplot(2,3,6)
    tSNEAbi2 = tsne(featureMatrixAbi2);
    gscatter(tSNEAbi2(:,1),tSNEAbi2(:,2),siteIdx, '','',10, 'off');
    title('SimS, SimE, & disHat(Z-scored)')

    [coeffJCom,scoreJCom,~] = pca(featureMatrixJCom);
    [coeffCom,scoreCom,~] = pca(featureMatrixCom);
    [coeffComI,scoreComI,~] = pca(featureMatrixComI);
    [coeffComZ,scoreComZ,~] = pca(featureMatrixComZ);
    [coeffXcorrZ,scoreXcorrZ,~] = pca(featureMatrixXcorrZ);
    [coeffAbi1,scoreAbi1,~] = pca(featureMatrixAbi1);
    [coeffAbi2,scoreAbi2,~] = pca(featureMatrixAbi2);

    figure(502);
    subplot(2,3,1)
    gscatter(scoreCom(:,1),scoreCom(:,2),siteIdx, '','',10, 'off');
    title('SimA, SimS, SimE, disHat, & NmAll (Raw)')
    subplot(2,3,2)
    gscatter(scoreComI(:,1),scoreComI(:,2),siteIdx, '','',10, 'off');
    title('SimA, SimS, SimE, disHat, & NmAll (Ranked)')
    subplot(2,3,3)
    gscatter(scoreComZ(:,1),scoreComZ(:,2),siteIdx, '','',10, 'off');
    title('SimA, SimS, SimE, disHat, & NmAll (Z-scored)')
    subplot(2,3,4)
    gscatter(scoreXcorrZ(:,1),scoreXcorrZ(:,2),siteIdx, '','',10, 'off');
    title('SimA, SimS, & SimE (Z-scored)')
    subplot(2,3,5)
    gscatter(scoreAbi1(:,1),scoreAbi1(:,2),siteIdx, '','',10, 'off');
    title('SimS, SimE, disHat, & NmAll (Z-scored)')
    subplot(2,3,6)
    gscatter(scoreAbi2(:,1),scoreAbi2(:,2),siteIdx, '','',10, 'off');
    title('SimS, SimE, & disHat(Z-scored)')

    while 0
        wanderlustRJCom=wanderlust(featureMatrixJCom);
        wanderlustRCom=wanderlust(featureMatrixCom);
        wanderlustRComI=wanderlust(featureMatrixComI);
        wanderlustRComZ=wanderlust(featureMatrixComZ);
        wanderlustRXcorrZ=wanderlust(featureMatrixXcorrZ);
        wanderlustRAbi1=wanderlust(featureMatrixAbi1);
        wanderlustRAbi2=wanderlust(featureMatrixAbi2);
    end

    wdrJCom = mean(wanderlustRJCom.T);
    wdrCom = mean(wanderlustRCom.T);
    wdrComI = mean(wanderlustRComI.T);
    wdrComZ = mean(wanderlustRComZ.T);
    wdrXcorrZ = mean(wanderlustRXcorrZ.T);
    wdrAbi1 = mean(wanderlustRAbi1.T);
    wdrAbi2 = mean(wanderlustRAbi2.T);

    [~,IwdrJCom] = sort(wdrJCom);
    [~,IwdrCom] = sort(wdrCom);
    [~,IwdrComI] = sort(wdrComI);
    [~,IwdrComZ] = sort(wdrComZ);
    [~,IwdrXcorrZ] = sort(wdrXcorrZ);
    [~,IwdrAbi1] = sort(wdrAbi1);
    [~,IwdrAbi2] = sort(wdrAbi2);

    pUpper = @(x) prctile(x,90);
    IwdrComUpper = splitapply(pUpper,IwdrCom,(siteIdx-1)/30+1);
    IwdrComIUpper = splitapply(pUpper,IwdrComI,(siteIdx-1)/30+1);
    IwdrComZUpper = splitapply(pUpper,IwdrComZ,(siteIdx-1)/30+1);
    IwdrAbi1Upper = splitapply(pUpper,IwdrAbi1,(siteIdx-1)/30+1);

    pLower = @(x) prctile(x,10);
    IwdrComLower = splitapply(pLower,IwdrCom,(siteIdx-1)/30+1);
    IwdrComILower = splitapply(pLower,IwdrComI,(siteIdx-1)/30+1);
    IwdrComZLower = splitapply(pLower,IwdrComZ,(siteIdx-1)/30+1);
    IwdrAbi1Lower = splitapply(pLower,IwdrAbi1,(siteIdx-1)/30+1);

    figure(7025)
    jbfill(unique(siteIdx),IwdrAbi1Upper,IwdrAbi1Lower,'g')
    hold on
    %jbfill(unique(siteIdx),IwdrComIUpper,IwdrComILower,'g');
    jbfill(unique(siteIdx),IwdrComZUpper,IwdrComZLower,'r');
    hold off
    title('Wanderlust (Z-scored)')
    ylabel('Rank')
    xlabel('Time (second)')
    xticks(time1)
    xticklabels(time1/30*0.2)
    legend('SimS, SimE, disHat, & NmAll', 'SimA, SimS, SimE, disHat, & NmAll');


    figure(7024)
    jbfill(unique(siteIdx),IwdrComUpper,IwdrComLower,'g')
    hold on
    %jbfill(unique(siteIdx),IwdrComIUpper,IwdrComILower,'g');
    jbfill(unique(siteIdx),IwdrComZUpper,IwdrComZLower,'r');
    hold off
    title('Wanderlust (SimA, SimS, SimE, disHat, & NmAll)')
    ylabel('Rank')
    xlabel('Time (second)')
    xticks(time1)
    xticklabels(time1/30*0.2)
    legend('Raw','Z-scored');

    figure(802024)
    %scatter(siteIdx,IwdrCom,'r')
    hold on
    %scatter(siteIdx,IwdrComI,'g');
    %scatter(siteIdx,IwdrComZ,'b');
    %scatter(siteIdx,IwdrXcorrZ,'b');
    scatter(siteIdx,IwdrAbi1,'r');
    %scatter(siteIdx,IwdrAbi2,'b');
    scatter(siteIdx,IdisHat,'b');
    hold off
    title('Wanderlust')
    ylabel('Arbitrary unit')
    xlabel('Time (second)')
    xticks(time1)
    xticklabels(time1/30*0.2)

    pcorJCom = corr(siteIdx',siteIdx(IwdrJCom)', 'Type', 'Spearman');
    pcorCom = corr(siteIdx',siteIdx(IwdrCom)', 'Type', 'Spearman');
    pcorComI = corr(siteIdx',siteIdx(IwdrComI)', 'Type', 'Spearman');
    pcorComZ = corr(siteIdx',siteIdx(IwdrComZ)', 'Type', 'Spearman');
    pcorDisHat = corr(siteIdx',siteIdx(IdisHat)', 'Type', 'Spearman');
    pcorXcorrZ = corr(siteIdx',siteIdx(IwdrXcorrZ)', 'Type', 'Spearman');
    pcorAbi1 = corr(siteIdx',siteIdx(IwdrAbi1)', 'Type', 'Spearman');
    pcorAbi2 = corr(siteIdx',siteIdx(IwdrAbi2)', 'Type', 'Spearman');
end
function y = findMinCor(x)
    global sites
    xnmrot1G1 = [];
    ynmrot1G1 = [];
    xnmrot2G1 = [];
    ynmrot2G1 = [];
    xnmrot1G2 = [];
    ynmrot1G2 = [];
    xnmrot2G2 = [];
    ynmrot2G2 = [];

    idxG1=find(x);
    idxG2=find(~x);
    for i = idxG1
        xnmrot1G1 = [xnmrot1G1; sites(i).evaluation.CME2CSide_yule2.locs1.xnmrot];
        ynmrot1G1 = [ynmrot1G1; sites(i).evaluation.CME2CSide_yule2.locs1.ynmrot];
        xnmrot2G1 = [xnmrot2G1; sites(i).evaluation.CME2CSide_yule2.locs2.xnmrot];
        ynmrot2G1 = [ynmrot2G1; sites(i).evaluation.CME2CSide_yule2.locs2.ynmrot];
    end

    for i = idxG2
        xnmrot1G2 = [xnmrot1G2; sites(i).evaluation.CME2CSide_yule2.locs1.xnmrot];
        ynmrot1G2 = [ynmrot1G2; sites(i).evaluation.CME2CSide_yule2.locs1.ynmrot];
        xnmrot2G2 = [xnmrot2G2; sites(i).evaluation.CME2CSide_yule2.locs2.xnmrot];
        ynmrot2G2 = [ynmrot2G2; sites(i).evaluation.CME2CSide_yule2.locs2.ynmrot];
    end
    G1 = getKernelMatrix(true, [500 500], [xnmrot1G1 ynmrot1G1; xnmrot2G1 ynmrot2G1]);
    G2 = getKernelMatrix(true, [500 500], [xnmrot1G2 ynmrot1G2; xnmrot2G2 ynmrot2G2]);
    
    D = xcorr2(G1,G2);
    y=log10(max(D(:)));
end

function [z, bw] = getKernelMatrix(centralized,matrixSize, sVrot, varargin) 
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

function groupDMap = calGruopRefDMap(idx, rangeV)
    global sites
    xnmrot1G = []; ynmrot1G = []; xnmrot2G = []; ynmrot2G = [];
    for i = rangeV
        xnmrot1G = [xnmrot1G; sites(idx(i)).evaluation.CME2CSide_yule2.locs1.xnmrot];
        ynmrot1G = [ynmrot1G; sites(idx(i)).evaluation.CME2CSide_yule2.locs1.ynmrot];
        xnmrot2G = [xnmrot2G; sites(idx(i)).evaluation.CME2CSide_yule2.locs2.xnmrot];
        ynmrot2G = [ynmrot2G; sites(idx(i)).evaluation.CME2CSide_yule2.locs2.ynmrot];
    end
    groupDMap = getKernelMatrix(true, [500 500], [xnmrot1G ynmrot1G; xnmrot2G ynmrot2G]); % Todo: the function to generate the image can be changed
end