%% Tempalte-free particle fusion of NPC using LocMoFit
% This script allows one to perform tempalte-free particle fusion using LocMoFit. Please make a copy before using this script and only work on the copy. Please define the parameters in the section of **Basic parameters** per task. After running through the analysis, you can find a the final average in the variable *finalAvg* in the _avg.mat file in ../analysis. finalAvg is a cell contains averages yeilded after all iterations. In the _avg.mat, the variable indBest tells you that finalAvg{indBest} is the best avarege, according to its mean log-likelihood.
%
% Output files:
%   * [anyIDYouLike]_preReg.mat: intermediate results of pre-registered pores.
%   * [anyIDYouLike]_LLMatrix.mat: intermediate results of the log-likelihood matrix.
%   * [anyIDYouLike]_initialTemp.mat: intermediate results of the initial template.
%   * [anyIDYouLike]_avg.mat: NPC averages.
%
% Last update:
%   20.12.2022
%
% Version:
%   v0.0.1
%
% Source:
%   template_free_particle_fusion.m
%
% What's new
%   The pre-registeration based on the dual-ring model is now disabled by default.

%% Basic parameters
% Note: the following parameters are optimized. Do not change them
% unless necessary

% source _sml.mat file related
par.sml_source = 'anyID';             % define an optional ID to represent the source _sml.mat file. This is not used in the following script.
par.sml_driveLoad = 'd';                     % the letter of the drive where the source data is.
par.sml_dirPath = ':\path\to\the\folder\where\the\file\is\saved\'; % the path to the folder.
par.sml_fileName = '_theFile_sml.mat';       % the name of the _sml.mat file.

% data save related
% locprecnm: localization presicion in the unit of nanometer.
par.save_wd = ':\template\analysis\';       % working directory without its drive letter
par.save_drive = 'd';                       % the letter of the drive where the working directory is.
par.save_name = 'anyIDYouLike';             % define the ID (used when saving) of the analysis here.

% parameters of the particle fusion
par.parFu_rng = 211023;                     % seed for the RNG.
par.parFu_path2Mod = "c:\whereTheSMAPFolderIs\SMAP\LocMoFit\models\NPCPointModel_flexible2.m"; % path to the dual-ring model.


%% Advanced parameters
% Note: the following parameters are optimized. Do not change them
% unless necessary

% parameters of the particle fusion
par.parFu_numOfInitSite = 50;               % a subset of particles used for the all-to-all registration.
par.parFu_maxIter = 100;                    % maximum number of iterations for the iterative registrations.
par.parFu_maxStall = 10;                    % maximum stall of improvement on the LL.
par.parFu_minImprov = 0.1;                  % maximum effective improvement of the LL.
par.parFu_method_initSite = 'sumLL';        % metric used for building the intial model. Can be either 'sumLL' or 'sumRank'.
par.parFu_useIsoLocprecnm = true;           % whether to use isotropic locprecnm (localization precision).
par.parFu_correctOrientation = false;       % correct the orientation the average NPCs everytime after fusion. 
par.parFu_preRegistration = false;           % determine the initial parameters of the registration based on the fit of the dual-ring model.

% settings for the procedure control
par.proc_autoFileLoading = true;            % whether to load the _sml.mat file automatically.
par.proc_continueOldAvg = true;             % whether to continue the analysis from an old averaging.
par.proc_arunSummarizeModFitNPC3D = false;  % whether to run SummarizeModFitNPC3D.


%% ---------Typically, users do not have to change the following code------------
%% Get started
clearvars -except par
clc
close all

%% Loading file
% Initiate SMAP and load an sml.mat file base on parameters par.sml_XXX. This
% creates an SMAP session as the object g, with the sml file loaded.
g = initParFu(par);

%% Initiate procedure control
% Detecting the previously saved analysis and deciding where to continue
% in the procedure. This creates a new field *proc_continueFrom* in the
% structure array *par*
par = initProcCtrl(par);

%% Correct tilt and tip.
disp('Correctting tilts and tips...')
[selectedSetInd,originalLocs, nUsedForAvg, preRegLocs, lPars_preReg, preRegFitPar_raw] = dataPrep(par, g);
save([par.save_drive par.save_wd par.save_name '_preReg.mat'], 'originalLocs', 'preRegLocs', 'lPars_preReg', 'preRegFitPar_raw')

%% All-to-all registration
disp('All-to-all registration starts...')
[LL, modelPars, startSet] = all2allInitialSet(par, originalLocs, preRegLocs, lPars_preReg);
save([par.save_drive par.save_wd par.save_name '_LLMatrix.mat'], 'startSet', 'selectedSetInd', 'modelPars', 'LL')

%% Determine the initial site
disp('Determine the initial site...')
[idxStart, rank_sum] = getInitSites(par, LL);

%% Setup of LocMoFit for generating the first data-driven template
disp('Building the first data-driven template...')
refLocs_ori = startSet{idxStart};
refLocs = fuseLocsNearby(refLocs_ori,2);
if par.parFu_preRegistration
    fitting_step2 = setUp_p2pFit(refLocs, 'eps', 5);
else
    fitting_step2 = setUp_p2pFit_free3Drot(refLocs, 'eps', 5);
end
fitting_step2.modelVerCascade = 1:3;

%% Show initial site
hCumulative = initGuiTabGrp('Cumulative fusion');
ax = initGuiTab(hCumulative, 'Site 1');
plotSite(ax, startSet{idxStart}, fitting_step2)
% plotSite(startSet{idxStart}, fitting_step2)

%% Build the initial template
[allNodes, all2Ref_Par] = buildInitTemplate(par, fitting_step2, refLocs_ori, refLocs, originalLocs, rank_sum, hCumulative, lPars_preReg);

%% Initializing the Iterative registration
disp('Interative registration starts...')
[refLocs, iter, iterStall, fitting_step2] = initIterReg(par, allNodes, rank_sum);

%% Running the iterative registration
[finalAvg, all2Avg_Par, all2Avg_LL, indBest] = iterReg(par, fitting_step2, refLocs, originalLocs, nUsedForAvg, iter, iterStall, lPars_preReg);

%% Load the final average to SMAP
loadAvg2SMAP(g, fitting_step1, fitting_step2, indBest, selectedSetInd,all2Avg_Par)
disp('Particle fusion done.')


%% Functions for each step
function g = initParFu(par)
% initParFu() initiates the particle fusion analysis.
%
% Usage:
%   g = initParFu(par)
%
% Args:
%   par: a struct of parameters of the analysis.
%
% Returns:
%   g: a GuiMainSMAP object.
% 

if par.proc_autoFileLoading
    % start a new SMAP session
    SMAP
    
    % set loading settings
    guiFile = g.children.guiFile;
    guiFile.guihandles.loadmodule.Value = 1;
    guiFile.guihandles.loadmodule.Callback{1}(guiFile.guihandles.loadmodule,[],guiFile)
    
    loader_1 = guiFile.children.loader_1;
    p = loader_1.getAllParameters;
    p.updateGuiPar = 1;
    loader_1.setGuiParameters(p);
    
    % load files
    rootPath = [par.sml_driveLoad par.sml_dirPath];
    lAdd = 0;
    g.children.guiFile.loadbutton_callback([],[],lAdd,rootPath,par.sml_fileName)
    
    % set locs filtering
    L1 = g.children.guiRender.children.Layer1;
    layerp=L1.getPar(L1.layerprefix);
    g.locData.setPar('layer1_selectedField', {'znm', layerp.znm_min, layerp.znm_max, 1, 1})
    g.locData.setPar('layer1_selectedField', {'frame', layerp.frame_min, layerp.frame_max, 1, 1})
    g.locData.setPar('layer1_selectedField', {'locprecnm', layerp.locprecnm_min, layerp.locprecnm_max, 1, 1})
    g.locData.setPar('layer1_selectedField', {'LLrel', layerp.LLrel_min, layerp.LLrel_max, 1, 1})
    layerp=L1.getPar(L1.layerprefix);
    L1.setPar(L1.layerprefix, layerp)
    L1.updateLayerField 
end
end

function par = initProcCtrl(par)
% initProcCtrl() initates the process control.
%
% Usage:
%   g = initProcCtrl(par)
%
% Args:
%   par: a struct of parameters of the analysis.
%
% Returns:
%   par: an updated struct of parameters of the analysis.
% 
if par.proc_continueOldAvg
    if exist([par.save_drive par.save_wd par.save_name '_avg.mat'],'file')==2
        par.proc_continueFrom = 3;
        load([par.save_drive par.save_wd par.save_name '_initialTemp.mat'])
        load([par.save_drive par.save_wd par.save_name '_LLMatrix.mat'])
        load([par.save_drive par.save_wd par.save_name '_preReg.mat'])
        load([par.save_drive par.save_wd par.save_name '_avg.mat'])
    elseif exist([par.save_drive par.save_wd par.save_name '_initialTemp.mat'],'file')==2
        par.proc_continueFrom = 2;
        load([par.save_drive par.save_wd par.save_name '_initialTemp.mat'])
        load([par.save_drive par.save_wd par.save_name '_LLMatrix.mat'])
        load([par.save_drive par.save_wd par.save_name '_preReg.mat'])
        disp('No old average found. Starting from the iteration 0.')
    elseif exist([par.save_drive par.save_wd par.save_name '_LLMatrix.mat'],'file')==2
        par.proc_continueFrom = 1;
        load([par.save_drive par.save_wd par.save_name '_LLMatrix.mat'])
        load([par.save_drive par.save_wd par.save_name '_preReg.mat'])
        disp('No initial template found. Starting from the initial template.')
    elseif exist([par.save_drive par.save_wd par.save_name '_preReg.mat'],'file')==2
        load([par.save_drive par.save_wd par.save_name '_preReg.mat'])
        par.proc_continueFrom = 0;
        disp('No pre-registration found. Starting from pre-registration.')
    else
        par.proc_continueFrom = -1;
        warning('Nothing to continue. Starting from scratch.')
    end
    variableNames = who;
    variableNames = variableNames(~strcmp(variableNames,'par'));
    for k=1:length(variableNames)
        assignin('base', variableNames{k}, eval(variableNames{k}))
    end
else
    par.proc_continueFrom = 0;
end
end
 
function [selectedSetInd,originalLocs, nUsedForAvg, preRegLocs, lPars_preReg, preRegFitPar_raw] = dataPrep(par, g)
% To-do: remove the QC part even
se = g.locData.SE;
% remove bad sites (one-ring or out-of-focus pores) detected in the
% previous analysis
if par.proc_arunSummarizeModFitNPC3D
    h_summary = openPlugin(g, {'ROIManager','Analyze','summarizeModFitNPC3D'});
    h_summary.run(h_summary.getAllParameters)
end
sites = se.sites;
list3Val = getFieldAsVector(sites,'annotation.list3.value');
lUse = getFieldAsVector(sites,'annotation.use');
lGood = list3Val==1&lUse;
nUsedForAvg = sum(lGood);

%% Randomly picked init sites
rng(par.parFu_rng)
selectedSets_ind = find(lGood);
selectedSetInd = datasample(selectedSets_ind, nUsedForAvg,'Replace',false); % randomly select pores

timerVal = tic;

if par.parFu_preRegistration
    preRegLocMoFitter = setUp_preRegFit(par);
else
    preRegLocMoFitter = [];
end

preRegLocs = [];
lPars_preReg = [];
preRegFitPar_raw = [];

if par.proc_continueFrom<0

    for k = length(selectedSetInd):-1:1
        % get locs for all sites
        indSite = selectedSetInd(k);
        pos = sites(indSite).pos;
        file = sites(indSite).info.filenumber;
        locs(k) = g.locData.getloc({'xnm','ynm','znm','locprecnm','channel','locprecznm','layer'},'grouping','grouped','position',[pos(1:2) 150 150],'filenumber',file, 'layer',1);
        locs(k).xnm = locs(k).xnm-pos(1);
        locs(k).ynm = locs(k).ynm-pos(2);
        % assign layer
        locs(k).layer = zeros(size(locs(k).xnm));
        locs(k).layer(locs(k).channel==0) = 1;

        if par.parFu_useIsoLocprecnm
            locs(k).locprecznm = locs(k).locprecnm;
        end
    end
    if ~license('test','Distrib_Computing_Toolbox')||length(selectedSetInd)<10
        for k = length(selectedSetInd):-1:1
            if par.parFu_preRegistration
                [preRegLocs{k}, lPars_preReg{k}, preRegFitPar_raw{k}] = preReg(preRegLocMoFitter,locs(k));
            else
                locs(k).xnm = locs(k).xnm - median(locs(k).xnm);
                locs(k).ynm = locs(k).ynm - median(locs(k).ynm);
                locs(k).znm = locs(k).znm - median(locs(k).znm);
            end

            originalLocs{k} = locs(k);
            %     selectedSet{k} = fitting_step1.locsHandler(locs, lParM1,1);

            %
            if toc(tic)>=5 % update every 15 sec
                timerVal = tic;
                disp(['Site ' num2str(k) ' done!'])
            end
        end
    else
        if par.parFu_preRegistration
            preRegLocs = cell(1,length(selectedSetInd));
            lPars_preReg = cell(1,length(selectedSetInd));
            preRegFitPar_raw = cell(1,length(selectedSetInd));
        end
        originalLocs = cell(1,length(selectedSetInd));
        parfor k = 1:length(selectedSetInd)
            if par.parFu_preRegistration
                [preRegLocs{k}, lPars_preReg{k}, preRegFitPar_raw{k}] = preReg(preRegLocMoFitter,locs(k));
            else
                locs(k).xnm = locs(k).xnm - median(locs(k).xnm);
                locs(k).ynm = locs(k).ynm - median(locs(k).ynm);
                locs(k).znm = locs(k).znm - median(locs(k).znm);
            end

            originalLocs{k} = locs(k);
            %     selectedSet{k} = fitting_step1.locsHandler(locs, lParM1,1);

            %
            disp(['Site ' num2str(k) ' done!'])
        end
    end
else
    originalLocs = evalin('base', 'originalLocs');
    preRegLocs = evalin('base', 'preRegLocs');
    lPars_preReg = evalin('base', 'lPars_preReg');
    preRegFitPar_raw = evalin('base', 'preRegFitPar_raw');
end
end

function [LL, modelPars, refSource] = all2allInitialSet(par, originalLocs, preRegLocs, lPars_preReg)
% function [LL, modelPars, fitInfo, startSet] = all2allInitialSet(par, selectedSet, fitting_step1, par_alignment, originalLocs)

%init
if isempty(preRegLocs)
    refSource = originalLocs(1:par.parFu_numOfInitSite);     % selectedSet was randomized.
    targetSource = refSource;
    lPars_preReg_startSet = [];
else
    refSource = preRegLocs(1:par.parFu_numOfInitSite);     % selectedSet was randomized.
    targetSource = originalLocs(1:par.parFu_numOfInitSite);
    lPars_preReg_startSet = lPars_preReg(1:par.parFu_numOfInitSite);     % selectedSet was randomized.
end
if par.proc_continueFrom < 1
    LL.iter1 = zeros(length(targetSource));      % log-likelihood.
    modelPars.iter1 = cell(length(targetSource));     % fitted parameters.
%     fitInfo.iter1 = cell(length(startSet));  % additional information.

    orderInMatrix = 1:length(targetSource);
    for l=1:length(targetSource)
        % l for the model
        locsRefOri = refSource{l};
        locsRef = locsRefOri;
        if l==1
            if par.parFu_preRegistration
                fitting_step2 = setUp_p2pFit(locsRef);
            else
                fitting_step2 = setUp_p2pFit_free3Drot(locsRef);
            end
        end
        
        ind = orderInMatrix(orderInMatrix~=l);% for sites
        % run particle to particle fitting
        [thisPars,thisLL] = p2pFit(fitting_step2, locsRef,targetSource, ind,'lPars_preReg', lPars_preReg_startSet, 'eps', 5);
        modelPars.iter1(l,ind) = thisPars;
%                 fitInfo.iter1{l,k} = fitting_step2.fitInfo;
        LL.iter1(l,ind) = thisLL;
        %             disp(['Site' num2str(k) ' to site' num2str(l) ' done!'])
        
        % for diagonal
        fitting_step2.fit(locsRef, 'controlLogLikelihood','overfitted');
        LL.iter1(l,l) = fitting_step2.fitInfo.LLOF;
    end
else
    LL = evalin('base', 'LL');
    modelPars = evalin('base', 'modelPars');
%     fitInfo = evalin('base', 'fitInfo');
end
end


function [idxStart, rank_sum] = getInitSites(par, LL)
switch par.parFu_method_initSite
    case 'sumRank'
        LL_iter1 = LL.iter1;
        LL_iter1_2 = LL_iter1;
        for k = 1:length(LL_iter1)
            LL_iter1(k,k) = -inf;
        end
        [~,idx] = sort(LL_iter1,1,'descend');
        r = 1:length(LL_iter1);
        rank = zeros(size(LL_iter1));
        for k = 1:length(LL_iter1)
            rank(idx(:,k),k) = r;
        end
        score = -sum(rank,2);% the higher the better
        [~,rank_sum] = sort(score,'descend');
        idxStart = rank_sum(1);
    case 'sumLL'
        LL_iter1 = LL.iter1;
        LL_iter1_2 = LL_iter1;
        for k = 1:length(LL_iter1)
            LL_iter1(k,k) = 0;
        end
        score = sum(LL_iter1,2);
        [~,rank_sum] = sort(-score);
        idxStart = rank_sum(1);
end
end

function [allNodes, all2Ref_Par] = buildInitTemplate(par, fitting_step2, refLocs_ori, firstRefLocs, originalLocs, rank_sum, hCumulative, lPars_preReg)
if par.proc_continueFrom < 2
    % init
    refLocs = firstRefLocs;
    allNodes = cell(50,1);
    all2Ref_Par = cell(50,1);

    for l=2:length(rank_sum)
        
        %     plotSite(allNodes{idxStart,idxOneAbove}, fitting_step2, false)
        
        idxTarget = rank_sum(l);
        
        if isempty(lPars_preReg)
            [modelPar, ~] = p2pFit(fitting_step2, refLocs,originalLocs, idxTarget);
        else
            [modelPar, ~] = p2pFit(fitting_step2, refLocs,originalLocs, idxTarget, 'lPars_preReg', lPars_preReg, 'eps',5);
        end
        
        locs = transByp2pFit(fitting_step2, originalLocs{idxTarget}, modelPar{1});
        refLocs_ori = fuseParticle(refLocs_ori, locs);
        allNodes{idxTarget} = refLocs_ori;
        refLocs = fuseLocsNearby(refLocs_ori,2);
        all2Ref_Par{idxTarget} = modelPar;
        
        ax = initGuiTab(hCumulative, ['Site ' num2str(l)]);
        plotSite(ax, allNodes{idxTarget}, fitting_step2, false)
        drawnow
        save([par.save_drive par.save_wd par.save_name '_initialTemp.mat'], 'allNodes','all2Ref_Par','rank_sum')
    end
else
    allNodes = evalin('base', 'allNodes');
    all2Ref_Par = evalin('base', 'all2Ref_Par');
    for l=2:length(rank_sum)
        idxTarget = rank_sum(l);
        ax = initGuiTab(hCumulative, ['Site ' num2str(l)]);
        plotSite(ax, allNodes{idxTarget}, fitting_step2, false)
        drawnow
    end
    
end
end


function [refLocs, iter, iterStall, fitting_step2] = initIterReg(par, allNodes, rank_sum)
if par.proc_continueFrom < 3
    refLocs = fuseLocsNearby(allNodes{rank_sum(end)},7);
    if par.parFu_preRegistration
        fitting_step2 = setUp_p2pFit(refLocs,'eps',5);
    else
        fitting_step2 = setUp_p2pFit_free3Drot(refLocs,'eps',5);
    end
    refLocs(2) = fuseLocsNearby(allNodes{rank_sum(end)},5);
    refLocs(3) = fuseLocsNearby(allNodes{rank_sum(end)},2);
    iter = 1;
    iterStall = 0;
else
    finalAvg = evalin('base', 'finalAvg');
    iter = evalin('base', 'iter');
    iterStall = evalin('base', 'iterStall');
    refLocs = fuseLocsNearby(finalAvg{end},7);
    if par.parFu_preRegistration
        fitting_step2 = setUp_p2pFit(refLocs,'eps',5);
    else
        fitting_step2 = setUp_p2pFit_free3Drot(refLocs,'eps',5);
    end
    refLocs(2) = fuseLocsNearby(finalAvg{end},5);
    refLocs(3) = fuseLocsNearby(finalAvg{end},2);
end
end


function [finalAvg, all2Avg_Par, all2Avg_LL, indBest] = iterReg(par, fitting_step2, refLocs, originalLocs, nUsedForAvg, iter, iterStall, lPars_preReg)
fig = figure('Name','Maximum log-likelihood');
ax = axes(fig);
if par.proc_continueFrom==3
    finalAvg = evalin('base', 'finalAvg');
    all2Avg_Par = evalin('base', 'all2Avg_Par');
    all2Avg_LL = evalin('base', 'all2Avg_LL');
    indBest = evalin('base', 'indBest');
else
    indBest = 0;
end
while iter <=par.parFu_maxIter && iterStall <= par.parFu_maxStall
    fitting_step2.modelVerCascade = 1:3;
    
    % fitting individual sites with the average generated in the last step
    if isempty(lPars_preReg)
        [all2Avg_Par{iter}, all2Avg_LL{iter}] = p2pFit(fitting_step2, refLocs,originalLocs, 1:nUsedForAvg, 'dispalyStatus', true);
    else
        [all2Avg_Par{iter}, all2Avg_LL{iter}] = p2pFit(fitting_step2, refLocs,originalLocs, 1:nUsedForAvg, 'dispalyStatus', true, 'lPars_preReg', lPars_preReg, 'eps',5);
    end
    % creating the new average by aligning all the sites to the previous
    % average followed by concatenating all localizations.
    finalAvg{iter} = createAvg(fitting_step2, originalLocs, all2Avg_Par{iter}, 1:nUsedForAvg);
    % determine convergence
    if iter > 1
        meanLLStep_all = cellfun(@mean,all2Avg_LL); % the first one is for the initial template so doesn't count
        meanLLStep = meanLLStep_all(2:end);
        [bestLL,indBest] = max(meanLLStep);
        lImproved = meanLLStep(end)-bestLL>par.parFu_minImprov;
    end
    
    if iter<=2||lImproved
        iterStall = 0;
    else
        iterStall = iterStall + 1;
    end
    
    if iter > 1
        % show LL for each iter
        plot(ax, 0:(iter-1), meanLLStep_all);
        xlabel(ax, 'Iteration');
        ylabel(ax, 'Mean maximum log-likelihood');

        drawnow
    end
    
    iter = iter+1;
    save([par.save_drive par.save_wd par.save_name '_avg.mat'], 'finalAvg', 'all2Avg_Par','iter','all2Avg_LL','iterStall','indBest')
    
    % Get ready for the next step
    if iter <=par.parFu_maxIter && iterStall <= par.parFu_maxStall
        refLocs = fuseLocsNearby(finalAvg{iter-1},7);
        refLocs(2) = fuseLocsNearby(finalAvg{iter-1},5);
        refLocs(3) = fuseLocsNearby(finalAvg{iter-1},2);
    end
end
end


function loadAvg2SMAP(g, fitting_step1, fitting_step2, avg2show, selectedSetInd,all2Avg_Par)
g.locData.loc.xnm_fusion = zeros(size(g.locData.loc.xnm));
g.locData.loc.ynm_fusion = zeros(size(g.locData.loc.xnm));
g.locData.loc.znm_fusion = zeros(size(g.locData.loc.xnm));
offSetFromZero = 1000;
se = g.locData.SE;
for k = length(selectedSetInd):-1:1
    indSite = selectedSetInd(k);
    pos = se.sites(indSite).pos;
    file = se.sites(indSite).info.filenumber;
    
    [locs, locsInd] = g.locData.getloc({'xnm','ynm','znm','locprecnm','channel','locprecznm'},'grouping','ungrouped','position',[pos(1:2) 150 150],'filenumber',file);
    if k == length(selectedSetInd)
        finalParInd = false(size(locsInd));
    end
    finalParInd = finalParInd|locsInd;
    
    locs.xnm = locs.xnm-pos(1);
    locs.ynm = locs.ynm-pos(2);
    
    fitting_step2.allParsArg = all2Avg_Par{avg2show}{k};
    locs = fitting_step2.locsHandler(locs,fitting_step2.exportPars(1,'lPar'),1);

    g.locData.loc.xnm_fusion(locsInd) = locs.xnm;
    g.locData.loc.ynm_fusion(locsInd) = locs.ynm;
    g.locData.loc.znm_fusion(locsInd) = locs.znm;
    
    disp(['Site ' num2str(k) ' done!'])
end
[locs, ~] = g.locData.getloc({'xnm_fusion','ynm_fusion','znm_fusion','locprecnm','channel','locprecznm'},'grouping','ungrouped');
locs = subsetStruct(locs, finalParInd);
locs = RenameField(locs,{'xnm_fusion','ynm_fusion','znm_fusion','channel'},{'xnm','ynm','znm','layer'});
locs.layer = locs.layer+1;
locs = correctOrientation(fitting_step1, locs);
g.locData.loc.xnm_fusion(finalParInd) = locs.xnm + offSetFromZero;
g.locData.loc.ynm_fusion(finalParInd) = locs.ynm + offSetFromZero;
g.locData.loc.znm_fusion(finalParInd) = locs.znm;
end
%%
%
function locs = transByp2pFit(fitter2, locs, modelPar)
% transformation after the first fit
% locprecnm should not be changed
fitter2.allParsArg = modelPar;
lPar = fitter2.exportPars(1,'lPar');
lPar.variation = 0;
locs = fitter2.locsHandler(locs,lPar,1);
end

function [thisPars,thisLL, locs] = p2pFit(fitting_step2, locsRef,startSet, ind,varargin)
% todo: for multi-layer
p = inputParser;
p.addParameter('dispalyStatus',false)
p.addParameter('lPars_preReg',[])
p.addParameter('eps',10)
p.parse(varargin{:});
p = p.Results;

fn = fieldnames(locsRef);

for k = 1:length(locsRef)
    for n = 1:length(fn)
        locsRefL1(k).(fn{n}) = locsRef(k).(fn{n})(locsRef(k).layer == 1);
    end
    for n = 1:length(fn)
        locsRefL2(k).(fn{n}) = locsRef(k).(fn{n})(locsRef(k).layer == 2);
    end
end
fitting_step2.model{1}.pseudoModel = locsRefL1;

if ~license('test','Distrib_Computing_Toolbox')||length(ind)<10
    % For parallel computation
    for k = 1:length(ind)
        ink = ind(k);
        locs = startSet{ink};
        fitting_step2.resetInit;
        if ~isempty(p.lPars_preReg)
            fitting_step2.setParArg('m1.lPar.x', 'value', p.lPars_preReg{ink}.x)
            fitting_step2.setParArg('m1.lPar.y', 'value', p.lPars_preReg{ink}.y)
            fitting_step2.setParArg('m1.lPar.z', 'value', p.lPars_preReg{ink}.z)
            fitting_step2.setParArg('m1.lPar.xrot', 'value', p.lPars_preReg{ink}.xrot)
            fitting_step2.setParArg('m1.lPar.yrot', 'value',p.lPars_preReg{ink}.yrot)
        end
        fitting_step2.setParArg('m1.lPar.variation', 'value', p.eps, 'lb', -inf, 'ub', inf, 'min', 0, 'max', p.eps);
        fitting_step2.fit(locs);
        fitting_step2.setParArg('m1.lPar.x', 'value', fitting_step2.getVariable('modelPars.m1.lPar.x'));
        fitting_step2.setParArg('m1.lPar.y', 'value', fitting_step2.getVariable('modelPars.m1.lPar.y'));
        fitting_step2.setParArg('m1.lPar.z', 'value', fitting_step2.getVariable('modelPars.m1.lPar.z'));
        %             fitting_step2.fitInfo;
        thisPars{k} = fitting_step2.allParsArg;
        thisLL(k) = fitting_step2.fitInfo.LLfit;
        if p.dispalyStatus
            disp(['Site' num2str(k) ' done!'])
        end
    end
else
    subStartSet = startSet(ind);
    dispalyStatus = p.dispalyStatus;
    lPreReg = ~isempty(p.lPars_preReg);
    if lPreReg
        lPars_preReg = p.lPars_preReg(ind);
    else
        lPars_preReg = cell(1,length(ind));
    end
    parfor k = 1:length(ind)
        locs{k} = subStartSet{k};
        fitting_step2_ = fitting_step2;
        fitting_step2_.resetInit;
        if lPreReg
            fitting_step2_.setParArg('m1.lPar.x', 'value',lPars_preReg{k}.x)
            fitting_step2_.setParArg('m1.lPar.y', 'value',lPars_preReg{k}.y)
            fitting_step2_.setParArg('m1.lPar.z', 'value',lPars_preReg{k}.z)
            fitting_step2_.setParArg('m1.lPar.xrot', 'value',lPars_preReg{k}.xrot)
            fitting_step2_.setParArg('m1.lPar.yrot', 'value',lPars_preReg{k}.yrot)
        end
        fitting_step2_.setParArg('m1.lPar.variation', 'value', p.eps, 'lb', -inf, 'ub', inf, 'min', 0, 'max', p.eps);
        fitting_step2_.fit(locs{k});
        fitting_step2_.setParArg('m1.lPar.x', 'value', fitting_step2_.getVariable('modelPars.m1.lPar.x'));
        fitting_step2_.setParArg('m1.lPar.y', 'value', fitting_step2_.getVariable('modelPars.m1.lPar.y'));
        fitting_step2_.setParArg('m1.lPar.z', 'value', fitting_step2_.getVariable('modelPars.m1.lPar.z'));
        %             fitting_step2.fitInfo;
        thisPars{k} = fitting_step2_.allParsArg;
        thisLL(k) = fitting_step2_.fitInfo.LLfit;
        
        if dispalyStatus
            disp(['Site' num2str(k) ' done!'])
        end
    end
end
end


function [preRegLocs, lPars_preReg, preRegFitPar_raw] = preReg(preRegLocMoFitter, locs)
    preRegLocMoFitter.fit(locs)
    preRegFitPar_raw = preRegLocMoFitter.allParsArg;
    lPars_preReg = preRegLocMoFitter.exportPars(1,'lPar');
    lPars = lPars_preReg;
    lPars.variation = 0;
    preRegLocs = preRegLocMoFitter.locsHandler(locs, lPars,1);
end

function particle = createAvg(fitter2, siteSet, fittedPar, siteInd)
for k = 1:length(fittedPar)
    fitter2.allParsArg = fittedPar{k};
    locs = siteSet{siteInd(k)};
    lPars = fitter2.exportPars(1,'lPar');
    lPars.variation = 0;
    locs = fitter2.locsHandler(locs,lPars,1);
    if k == 1
        particle = locs;
    else
        particle = fuseParticle(particle, locs);
    end
end

% plotSite(particle, fitter2, false)
% drawnow;
% plotSite(particle, fitter2, true)
% drawnow;
end

function particle = correctOrientation(fitter, particle)
fitter.fit(particle);
lParM1 = fitter.exportPars(1,'lPar');
lParM2 = fitter.exportPars(2,'lPar');
lParM1.x = lParM1.x+fitter.rel(lParM2.z,3,1);
lParM1.y = lParM1.y+fitter.rel(lParM2.z,3,2);
lParM1.z = lParM1.z+fitter.rel(lParM2.z,3,3)/2;
particle = fitter.locsHandler(particle, lParM1,1);
end

function plotSite(varargin)
if isa(varargin{1}, 'matlab.graphics.axis.Axes')
    ax = varargin{1};
    particle = varargin{2};
    fitter = varargin{3};
else
    particle = varargin{1};
    fitter = varargin{2};
end
% plot when succeeded
% locsRefOri = particle;
locsRef = particle;


fn = fieldnames(locsRef);
for n = 1:length(fn)
    locsRefL1.(fn{n}) = locsRef.(fn{n})(locsRef.layer == 1);
end

fitter.model{1}.pseudoModel = locsRefL1;

if ~exist('ax','var')
    figure; tiledlayout(1,2);
    nexttile
    imagesc(squeeze(sum(fitter.model{1}.getImage([],'useLocprecnm',false, 'sigma', 5),3)))
    axis equal
    nexttile
    imagesc(squeeze(sum(fitter.model{1}.getImage([],'useLocprecnm',false, 'sigma', 5),1))')
    axis equal
else
    hTab = ax.Parent;
    clear('ax');
    t = tiledlayout(hTab, 1,2);
    hView = nexttile(t);
    imagesc(hView, squeeze(sum(fitter.model{1}.getImage([],'useLocprecnm',false, 'sigma', 5),3)))
    axis equal
    hView = nexttile(t);
    imagesc(hView, squeeze(sum(fitter.model{1}.getImage([],'useLocprecnm',false, 'sigma', 5),1))')
    axis equal
end
end

function [ID_target, ID_oneAbove] = findTargets(Z, ID)
% findTargets finds the targets of the site with the input ID.
% Arg:
%
nLeafNodes = size(Z,1)+1;
% find the current row
[row,col] = find(Z(:,1:2)==ID);
ID_samePair = Z(row, 3-col);
% find one node above
ID_oneAbove = nLeafNodes+row;
ID_target = [];
if ID_samePair>nLeafNodes
    descendant_samePair = Z(ID_samePair-nLeafNodes,1:2);
    ID_target = [ID_target descendant_samePair(descendant_samePair<=nLeafNodes)];
    haveDescendant = descendant_samePair(descendant_samePair>nLeafNodes);
    nHaveDescendant = length(haveDescendant);
    while nHaveDescendant>0
        haveDescendant_2 = [];
        for k = 1:nHaveDescendant
            descendant = Z(haveDescendant(k)-nLeafNodes,1:2);
            ID_target = [ID_target descendant(descendant<=nLeafNodes)];
            haveDescendant_2 = [haveDescendant_2 descendant(descendant>nLeafNodes)];
        end
        haveDescendant = haveDescendant_2;
        nHaveDescendant = sum(haveDescendant>nLeafNodes);
    end
else
    ID_target = ID_samePair;
end
end

function fitting_step2 = setUp_p2pFit(locsRef, varargin)
p = inputParser;
p.addParameter('eps',10);
p.parse(varargin{:});
p = p.Results;

mod = locsModel(locsRef, 'layer', 1);
mod.dimension = 3;
mod.layer = 1;

fitting_step2 = LocMoFit('SolverName', 'fminsearchbnd','SolverOptions',{'Display','off'});
fitting_step2.objFunType = 'likelihood';
fitting_step2.dataDim = 3;
fitting_step2.addModel(mod);
fitting_step2.sigmaCascade = [1 1 1; 0 0 0];
fitting_step2.model{1}.sigmaFactor = [1 0];

fitting_step2.roiSize = 300;
fitting_step2.setParArg('m1.lPar.xrot','fix',false,'lb',-25,'ub',25)
fitting_step2.setParArg('m1.lPar.yrot','fix',false,'lb',-25,'ub',25)
fitting_step2.setParArg('m1.lPar.variation','fix',false,'value',p.eps,'lb',-inf,'ub',inf,'min',0,'max',p.eps)
fitting_step2.setParArg('m1.lPar.zrot','fix',false,'lb',-180,'ub',180)
fitting_step2.setParArg('m1.lPar.x','lb',-25,'ub',+25)
fitting_step2.setParArg('m1.lPar.y','lb',-25,'ub',+25)
fitting_step2.setParArg('m1.lPar.z','lb',-25,'ub',+25,'value',0)
fitting_step2.setParArg('m1.lPar.weight','fix',true,'value',1)
fitting_step2.setParArg('m91.offset.weight','fix',false, 'value', 0.1,'lb',-inf,'ub',inf,'min',1e-3,'max',0.999)

fitting_step2.converter(fitting_step2, 'rand(1)*360', 'm1.lPar.zrot')

fitting_step2.advanceSetting.gaussDistMode.value = 'fast';
end

function fitting_step2 = setUp_p2pFit_free3Drot(locsRef, varargin)
p = inputParser;
p.addParameter('eps',10);
p.parse(varargin{:});
p = p.Results;

mod = locsModel(locsRef, 'layer', 1);
mod.dimension = 3;
mod.layer = 1;

fitting_step2 = LocMoFit('SolverName', 'fminsearchbnd','SolverOptions',{'Display','off'});
fitting_step2.objFunType = 'likelihood';
fitting_step2.dataDim = 3;
fitting_step2.addModel(mod);
fitting_step2.sigmaCascade = [1 1 1; 0 0 0];
fitting_step2.model{1}.sigmaFactor = [1 0];

fitting_step2.roiSize = 300;
fitting_step2.setParArg('m1.lPar.xrot','fix',false,'lb',-inf,'ub',inf,'min',-inf,'max',inf)
fitting_step2.setParArg('m1.lPar.yrot','fix',false,'lb',-inf,'ub',inf,'min',-inf,'max',inf)
fitting_step2.setParArg('m1.lPar.variation','fix',false,'value',p.eps,'lb',-inf,'ub',inf,'min',0,'max',p.eps)
fitting_step2.setParArg('m1.lPar.zrot','fix',false,'lb',-inf,'ub',inf,'min',-inf,'max',inf)
fitting_step2.setParArg('m1.lPar.x','lb',-25,'ub',+25)
fitting_step2.setParArg('m1.lPar.y','lb',-25,'ub',+25)
fitting_step2.setParArg('m1.lPar.z','lb',-25,'ub',+25,'value',0)
fitting_step2.setParArg('m1.lPar.weight','fix',true,'value',1)
fitting_step2.setParArg('m91.offset.weight','fix',false, 'value', 0.1,'lb',-inf,'ub',inf,'min',1e-3,'max',0.999)

fitting_step2.converter(fitting_step2, 'rand(1)*360', 'm1.lPar.zrot')

fitting_step2.advanceSetting.gaussDistMode.value = 'fast';
end

function preRegLocMoFitter = setUp_preRegFit(par)
mod = functionModel(par.parFu_path2Mod);
mod.dimension = 3;
mod.layer = 1;
mod.modelObj.setInternalSettings('copyPerCorner', 1)
mod.modelObj.setInternalSettings('cornerNum', 100)

preRegLocMoFitter = LocMoFit('SolverName', 'fminsearchbnd','SolverOptions',{'Display','off'});
preRegLocMoFitter.objFunType = 'likelihood';
preRegLocMoFitter.dataDim = 3;
preRegLocMoFitter.addModel(mod);

preRegLocMoFitter.roiSize = 250;

preRegLocMoFitter.setParArg('m1.lPar.xrot','fix',false,'lb',-25,'ub',25)
preRegLocMoFitter.setParArg('m1.lPar.yrot','fix',false,'lb',-25,'ub',25)
preRegLocMoFitter.setParArg('m1.lPar.variation','fix',false,'value',30,'lb',-inf,'ub',inf,'min',0,'max',30)
preRegLocMoFitter.setParArg('m1.lPar.zrot','fix',true,'value',0)
preRegLocMoFitter.setParArg('m1.lPar.x','lb',-30,'ub',+30)
preRegLocMoFitter.setParArg('m1.lPar.y','lb',-30,'ub',+30)
preRegLocMoFitter.setParArg('m1.lPar.z','lb',-30,'ub',+30)
preRegLocMoFitter.setParArg('m1.lPar.weight','fix',true,'value',1)
preRegLocMoFitter.setParArg('m91.offset.weight','fix',false, 'value', 0.1,'lb',-inf,'ub',inf,'min',1e-3,'max',0.999)
preRegLocMoFitter.setParArg('m1.mPar.ringDistance','lb',-inf,'ub',inf, 'value', 40,'min',0,'max',80)
preRegLocMoFitter.setParArg('m1.mPar.radius', 'value', 40,'lb',-inf,'ub',inf,'min',30,'max',70)

preRegLocMoFitter.converter(preRegLocMoFitter, 'median(locs.xnm)', 'm1.lPar.x')
preRegLocMoFitter.converter(preRegLocMoFitter, 'median(locs.ynm)', 'm1.lPar.y')
preRegLocMoFitter.converter(preRegLocMoFitter, 'median(locs.znm)', 'm1.lPar.z')
end