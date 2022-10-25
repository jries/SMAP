function out=run_LocMoFitGUI(obj, inp, varargin)
% This is the 'run' function for the LocMoFit GUI
% To-do:
%   To separate the GUI and LocMoFit layers.
%
eval_  = obj.locData.SE.processors.eval;
loadedPlugins = fieldnames(eval_.children);
allLocMoFitGUI = loadedPlugins(startsWith(loadedPlugins, 'LocMoFit'));

anyGUI_previewMode = 0;
for k = 1:length(allLocMoFitGUI)
    anyGUI_previewMode = anyGUI_previewMode + sum(eval_.children.(allLocMoFitGUI{k}).lPreview);
end

forceDisplay_currentGUI = obj.lPreview;
if anyGUI_previewMode
    % any one of steps is in the preview mode then set this
    % step also to the review mode
    varargin =  {'onlySetUp',true, 'forceDisplay',true, 'keepParsVal',true};
end

if strcmp(obj.fitter.getAdvanceSetting('siteSpecificMode'), 'on')&&~anyGUI_previewMode
    if isfield(obj.site.evaluation.(obj.name), 'written')
        obj.readPar
        inp = obj.getGuiParameters;
    else
        %                     do nothing
    end
end

if any(forceDisplay_currentGUI)
    % If this step is origianlly in the preview mode, throw
    % warning.
    warning off backtrace
    warning on verbose
    msg = ['Now is in the preview mode so the fit is not executed. To exit the preview mode, uncheck the preview box(es) of model(s) ' num2str(find(forceDisplay_currentGUI)) ' in ' obj.name];
    warning(msg, 'verbose')
    warning on backtrace
    warning off verbose
end

p = inputParser;
addParameter(p,'onlySetUp',false, @islogical);% not used at the moment
addParameter(p,'forceDisplay',false, @islogical); % only disply not fit
addParameter(p,'keepParsVal',false, @islogical); % keep the original pars
parse(p,varargin{:});
results = p.Results;
obj.fitter.linkedGUI = obj;
forceDisplay = results.forceDisplay;
keepParsVal = results.keepParsVal;

try
    %% RUNSMLMMODELFITGUI Actions to request from the plugin SMLMModelFitGUI
    % This is the 'run' function for the SMLMModelFitGUI.
    % Import the kimograph
    %             site = obj(1).locData

    %% Setup the models
    if ~inp.noFit||isempty(obj.fitter)||forceDisplay
        obj = setupModel(obj, inp);
    end

    fitter = obj.fitter;

    %% When other GUI are in the preview mode but not this one, do nothing but pass on the old output
    lPassOnDataOnly = ~any(forceDisplay_currentGUI) & anyGUI_previewMode;

    fetchPar(fitter, obj);

    %% Deal with locs
    [fieldQ, locs] = processLocs(obj);

    %% assign representative Sigma to each model
    % for m = 1:fitter.numOfModel
    %     thisLayer = fitter.model{m}.layer;
    %     lThisLayer = locs.layer==thisLayer;
    %     variation = fitter.getVariable(['pars.m' num2str(m) '.lPar.variation']);
    %     fitter.model{m}.locsMedSig = median(locs.locprecnm(lThisLayer).^2+variation.^2);
    % end
    %% Perform the fitting and export related
    if ~(inp.noFit||forceDisplay)
        %     fitter.resetInit;
        % also gets controlLogLikelihood
        fitter.fit(locs);
        % save the fit results
        out.allParsArg = fitter.allParsArg;
        out.parsInit = fitter.parsInit;
        out.fitInfo = fitter.fitInfo;
        out.externalInfo = fitter.externalInfo;
        obj.fitter = fitter;
        %         ax = fitter.plot(locs,'plotType','point');
    else
        % get info from the previous fit
        if ~forceDisplay
            % reuse fit results
            fitter.allParsArg = obj.locData.SE.currentsite.evaluation.(obj.name).allParsArg;
            fitter.parsInit = obj.locData.SE.currentsite.evaluation.(obj.name).parsInit;
            fitter.fitInfo = obj.locData.SE.currentsite.evaluation.(obj.name).fitInfo;
            if isfield(obj.locData.SE.currentsite.evaluation.(obj.name), 'externalInfo')
                fitter.externalInfo = obj.locData.SE.currentsite.evaluation.(obj.name).externalInfo;
            end
            fitter.fit(locs,'skipFit', true);
            out.fitInfo = fitter.fitInfo;
            out.allParsArg = fitter.allParsArg;
        else
            % preview
            tempAllParsArg = fitter.allParsArg;
            fitter.resetInit;
            % do conversion for preview
            fitter.convertNow(locs);
            if keepParsVal
                out.allParsArg = tempAllParsArg;
            end
        end
    end

    if isfield(obj.site.evaluation, obj.name)&&isfield(obj.site.evaluation.(obj.name), 'written')
        out.written = obj.site.evaluation.(obj.name).written;
    end
    %% Alignment (for averaging)
    if inp.useAlignment
        [obj, fitter] = registerLocs(obj, fieldQ, fitter);
    end
    %%
    %

    %% Display the result
    if inp.noFit||forceDisplay
        if keepParsVal
            tempAllParsArg = fitter.allParsArg; % this is to temporarily save the original allParsArg.
        end

        out.parsInit = fitter.parsInit;
        fitter.getDerivedPars;
        if ~strcmp(fitter.model{1}.modelType, 'image')
            fitter.deriveSigma(locs);
        end
        %% Get compensationFactor
        % This factor compensate the number of locs between different channels
        %                     compensationFactor = zeros(size(locs.layer))';
        %                     for k = 1:fitter.numOfLayer
        %                         compensationFactor(locs.layer==k) = fitter.compensationFactor(k);
        %                     end
        %                     compensationFactor(~ismember(locs.layer, fitter.allModelLayer)) = [];
        [~,~,parBestFit] = fitter.prepFit;
        %     fitter.fitInfo.ELL = fitter.getELL(parBestFit', compensationFactor, 2);

        out.fitInfo = fitter.fitInfo;
        out.externalInfo = fitter.externalInfo;
    end
    if ~lPassOnDataOnly
        if obj.getPar('se_display')||forceDisplay
            if fitter.model{1}.dimension == 2
                vis = obj.setoutput('Plot');
                fitter.showFitResult(vis, locs);
            else
                if ~strcmp(fitter.model{1}.modelType,'image')
                    viz1 = obj.setoutput('FreeRot');
                    viz2 = obj.setoutput('FixRot');
                    viz1Parent = viz1.Parent;
                    delete(viz1Parent.Children)
                    viz1 = axes(viz1Parent);

                    viz2Parent = viz2.Parent;
                    if isa(viz2Parent, 'matlab.graphics.layout.TiledChartLayout')
                        viz2Parent = viz2Parent.Parent;
                    end
                    delete(viz2Parent.Children)
                    for k = length(obj.locData.layer):-1:1
                        layerPar = obj.locData.SE.getLayerParameters(k,{'lut'});
                        allLut{k} = layerPar.lut.selection;
                        if k<=fitter.numOfModel
                            fitter.model{k}.displayLut = allLut{k};
                        end
                    end

                    fitter.plotFreeRot(viz1,locs,'lutLocs',allLut,'sigma',fitter.model{1}.sigma,'pixelSize',fitter.model{1}.pixelSize);
                    fitter.plotFixRot(viz2Parent,locs,'lutLocs',allLut,'sigma',fitter.model{1}.sigma,'pixelSize',fitter.model{1}.pixelSize);
                else
                    ax = fitter.plot(locs,'plotType','image','pixelSize',5);
                    ax2 = fitter.plot(locs,'plotType','image','pixelSize',5,'projection','yz');
                    ax3 = fitter.plot(locs,'plotType','image','pixelSize',5,'projection','xz');
                    vis = obj.setoutput('XY');
                    vis2 = obj.setoutput('YZ');
                    vis3 = obj.setoutput('XZ');
                    hFig1= ax.Parent;
                    hFig2= ax2.Parent;
                    hFig3= ax3.Parent;
                    ax.Parent = vis.Parent;
                    ax2.Parent = vis2.Parent;
                    ax3.Parent = vis3.Parent;
                    delete(vis)
                    delete(vis2)
                    delete(vis3)
                    close(hFig1)
                    close(hFig2)
                    close(hFig3)
                end
            end

            % Create the tab for fitted parameters
            axPar = obj.setoutput('Fitted_Par');

            % Clean the tab (remove existed uitable and so)
            delete(findobj(axPar.Parent.Children, '-not','type', 'axes'))
            uiFittedPar(axPar, fitter);
            axPar.Visible = 'Off';

            if keepParsVal
                fitter.allParsArg = tempAllParsArg; % since the display is only for this site, the original allParsArg should be put back.
            end
        end
    else
        for k = fitter.numOfModel:-1:1
            inp.(['pickSite_' num2str(k)]) = 0;
        end
        obj.setGuiParameters(inp)
    end
    out.fitInfo.guiInfo = 'Normal.';

    if strcmp(obj.fitter.getAdvanceSetting('siteSpecificMode'), 'on')&&~anyGUI_previewMode
        resetPar(obj)   
    end
catch ME
    if isfield(obj.site.evaluation, obj.name)
        out = obj.site.evaluation.(obj.name);
    end
    warning(['Model fitter did not run through. Site ' num2str(obj.site.ID) ' encountered some issues.'])
    disp(getReport(ME, 'extended', 'hyperlinks', 'on'))
    if results.forceDisplay
        out.fitInfo.guiInfo = 'Plot failed.';
    else
        out.fitInfo.guiInfo = 'Fit or plot failed.';
    end
    if strcmp(obj.fitter.getAdvanceSetting('siteSpecificMode'), 'on')&&~anyGUI_previewMode
        obj.resetPar
    end
end
end

function obj = setupModel(obj, inp)
%% Come to this part only when there is no the SMLMModelFit obj or when the user asks for.
fitter = obj.fitter;
if obj.getPar(['lOptimizerSettingEdited_' obj.name])
    % create the fitter obj.
    fitter.solver.SolverName = inp.optimizer.selection;      % set the solver
    optimizerparData = obj.guihandles.optimizerpar.Data';
    % add the settings to the solver
    optimizerparData = cellfun(@(x) iif(strcmp(x,'true'),@()true,~strcmp(x,'true'),@()x), optimizerparData,'UniformOutput',false);
    fitter.solver.SolverOptions={optimizerparData{:}};
    obj.setPar(['lOptimizerSettingEdited_' obj.name],false);
end
%         Deal with the model settings
if obj.getPar(['logicalParTableEdited_' obj.name])
    %% Set the model settings
    for k=obj.numMod:-1:1
        kStr = num2str(k);
        %                             fitter.model{k}.weight = inp.(['weight_fit_' kStr]); % kStr corresponds to the model order
        %             fitter.model{k}.layer = inp.(['layer_' kStr]); % this now
        %             will be updated right after editing
        if ~strcmp(fitter.model{k}.modelType,'image')
            fitter.model{k}.sigma = inp.(['sigma_fit_' kStr]);
            fitter.model{k}.sigmaFactor = inp.(['sigmaFactor_fit_' kStr]);
        end
        fitter.model{k}.pixelSize = inp.(['pixelsizefit_' kStr]);
        htab=obj.guihandles.(['partable_' kStr]); %handle of table
    end
    obj.setPar(['logicalParTableEdited_' obj.name], false);
end
fitter.roiSize = obj.P.par.se_siteroi.content;

% Converter
hConvert=obj.guihandles.anchorConvert; %handle of table
convertRules = hConvert.Data;
fitter.converterRules.target = {};
fitter.converterRules.rule = {};
for k = 1:size(convertRules,1)
    if strcmp(convertRules(k,1),'this step')
        sourceObj = [];
    else
        nameAllLoadedEval = obj.locData.SE.processors.eval.guihandles.modules.Data(:,2);
        indSourceObj = find(ismember(nameAllLoadedEval,convertRules(k,1)));
        sourceObj = obj.locData.SE.processors.eval.processors{indSourceObj}.fitter;
    end
    
    if length(convertRules{k,4})==0
        target = convertRules{k,3};
    else
        target = ['usr_' convertRules{k,4}];
    end
    
    fitter.converter(sourceObj,convertRules{k,2}, target);
end
obj.fitter = fitter;
end

function fetchPar(fitter, obj)
% Parameter arguments
for k = 1:fitter.numOfModel
    allParsArg = obj.guihandles.(['partable_' num2str(k)]).Data;
    allParsArg_size = size(allParsArg);
    for l = 1:allParsArg_size(1)
        parID = ['m' num2str(k) '.' strtrim(allParsArg{l,6}) '.' strtrim(allParsArg{l,1})];
        fitter.setParArg(parID, 'value', allParsArg{l,2}, 'fix', allParsArg{l,3},'lb', allParsArg{l,4},'ub', allParsArg{l,5}, 'min', allParsArg{l,7}, 'max', allParsArg{l,8}, 'label', allParsArg{l,9})
    end
end

% Offset arguments
offSetArg = obj.guihandles.layerSetting.Data;
offSetArg_size = size(offSetArg);
for l = 1:offSetArg_size(1)
    layerID = offSetArg{1,1}(end);
    BGForm = fitter.getAdvanceSetting(['m9' layerID '_background']);
    parID = ['m9' layerID '.offset.' BGForm];
    fitter.setParArg(parID, 'value', offSetArg{l,2}, 'fix', offSetArg{l,3},'lb', offSetArg{l,4},'ub', offSetArg{l,5}, 'min', offSetArg{l,6}, 'max', offSetArg{l,7}, 'label', '')
end
end

function [fieldQ, locs] = processLocs(obj)
layerson = find(obj.getPar('sr_layerson'));             % check the layers used
visitFlag = false;                                      % a flag for the conditional loop for the first round
if obj.fitter.dataDim == 3
    fieldQ = {'locprecnm','locprecznm','PSFxnm','xnm','znm','ynm','frame','xnmrot','ynmrot'};    % fields will be used
else
    fieldQ = {'locprecnm','PSFxnm','xnm','ynm','frame','xnmrot','ynmrot'};    % fields will be used
end
fieldQNLayer = [fieldQ 'layer'];
for k = layerson                                        % go through layers, collect all filtered locs
    if ~visitFlag
        locs=obj.getLocs(fieldQ,'layer',k,'size','freeroi');
        locs.layer = ones(size(locs.(fieldQ{1})))*k;
        fieldExact = fieldnames(locs);
        lRm = ~ismember(fieldExact, fieldQNLayer);
        locs = rmfield(locs,fieldExact(lRm));
    else
        locsNewLayer = obj.getLocs(fieldQ,'layer',k,'size','freeroi');
        for l = length(fieldQNLayer):-1:1
            if l == length(fieldQNLayer)
                locsNewLayer.layer = ones(size(locsNewLayer.(fieldQNLayer{l-1})))*k;
            end
            if isfield(locs, fieldQNLayer{l})
                locs.(fieldQNLayer{l})=[locs.(fieldQNLayer{l});locsNewLayer.(fieldQNLayer{l})];
            end
        end
    end
    fldElement = structfun(@(fld) length(fld),locs);
    fnLocs = fieldnames(locs);
    indFn2Rm = find(fldElement==0);
    locsOri = locs;
    locs = [];
    for l = length(fnLocs):-1:1
        if ~ismember(l,indFn2Rm)
            locs.(fnLocs{l}) = locsOri.(fnLocs{l});
        end
    end
    locs.xnm = locs.xnmrot;
    locs.ynm = locs.ynmrot;
    visitFlag = true;
    
end
end

function [obj, fitter] = registerLocs(obj, fieldQ, fitter)
% check the locs have been aligned or not
% if not, create empty fields for the new info after fitting
if ~isfield(obj.locData.loc, ['xnmaligned_' obj.name])
    obj.locData.loc.(['xnmaligned_' obj.name]) = single(zeros(size(obj.locData.loc.xnm)));
    obj.locData.loc.(['ynmaligned_' obj.name]) = single(zeros(size(obj.locData.loc.xnm)));
    obj.locData.loc.(['siteID_' obj.name]) = single(zeros(size(obj.locData.loc.xnm)));
    obj.locData.loc.(['rank_' obj.name]) = single(zeros(size(obj.locData.loc.xnm)));
    obj.locData.loc.class = single(zeros(size(obj.locData.loc.xnm)));
    if isfield(obj.locData.loc, 'znm')
        obj.locData.loc.(['znmaligned_' obj.name]) = single(zeros(size(obj.locData.loc.xnm)));
    end
end

% check a specific locs should be re-asigned or not based on the site ID
% if the site ID matches (saved one vs current one), then reset them.
lReset = obj.locData.loc.(['siteID_' obj.name])==obj.site.ID;
obj.locData.loc.(['siteID_' obj.name])(lReset) = 0;
obj.locData.loc.(['rank_' obj.name])(lReset) = 0;
obj.locData.loc.(['xnmaligned_' obj.name])(lReset) = 0;
obj.locData.loc.(['ynmaligned_' obj.name])(lReset) = 0;
obj.locData.loc.class(lReset) = 0;
if isfield(obj.locData.loc, 'znm')
    obj.locData.loc.(['znmaligned_' obj.name])(lReset) = 0;
end

% get locs
[newlocs, indNewLocs]=obj.getLocs(fieldQ,'size',obj.getPar('se_siteroi')/2,'group','ungrouped');
newlocs.xnm = newlocs.xnmrot;
newlocs.ynm = newlocs.ynmrot;

%% alignment
% convert
oldRules = fitter.converterRules;
% reset the rules
fitter.converterRules = [];
fitter.converterUserDefined = [];
alignSettings = obj.alignSettings;
for k = 1:size(alignSettings,1)
    if isempty(alignSettings{k,3})
        fitter.converter([],alignSettings{k,1},alignSettings{k,2})
    else
        fitter.converter([],alignSettings{k,1},['usr_' alignSettings{k,3}])
    end
end
fitter.convertNow
newlocs = fitter.locsRegister(newlocs, fitter.exportPars(1,'lPar'), 1);
fitter.converterRules = oldRules;

% save to the fields
idxNewLocs = find(indNewLocs);
lWithoutVal = obj.locData.loc.(['xnmaligned_' obj.name])(indNewLocs) == 0;
% postFitTrans do the transformation after corrected by the model
obj.locData.loc.(['xnmaligned_' obj.name])(idxNewLocs(lWithoutVal)) = newlocs.xnm(lWithoutVal)+500;
obj.locData.loc.(['ynmaligned_' obj.name])(idxNewLocs(lWithoutVal)) = newlocs.ynm(lWithoutVal)+500;
if isfield(newlocs,'znm')
    obj.locData.loc.(['znmaligned_' obj.name])(idxNewLocs(lWithoutVal)) = newlocs.znm(lWithoutVal);
end
obj.locData.loc.(['siteID_' obj.name])(idxNewLocs(lWithoutVal)) = obj.site.ID;
obj.locData.loc.class(idxNewLocs(lWithoutVal)) = obj.site.ID;
obj.locData.loc.(['rank_' obj.name])(idxNewLocs(lWithoutVal)) = obj.site.indList;

% if one loc appears in more than one site, assign the spatial info
% based on the distance to origin.
if sum(~lWithoutVal)>0
    idxWithVal = idxNewLocs(~lWithoutVal);
    if isfield(newlocs,'znm')
        distOriSq = obj.locData.loc.(['xnmaligned_' obj.name])(idxWithVal).^2+obj.locData.loc.(['ynmaligned_' obj.name])(idxWithVal).^2+obj.locData.loc.(['znmaligned_' obj.name])(idxWithVal).^2;
        distNewSq = newlocs.xnm(~lWithoutVal).^2+newlocs.ynm(~lWithoutVal).^2+newlocs.znm(~lWithoutVal).^2;
    else
        distOriSq = obj.locData.loc.(['xnmaligned_' obj.name])(idxWithVal).^2+obj.locData.loc.(['ynmaligned_' obj.name])(idxWithVal).^2;
        distNewSq = newlocs.xnm(~lWithoutVal).^2+newlocs.ynm(~lWithoutVal).^2;
    end
    idxSubWithVal = find(~lWithoutVal);
    lReplace = distNewSq>distOriSq;
    obj.locData.loc.(['xnmaligned_' obj.name])(idxWithVal(lReplace)) = newlocs.xnm(idxSubWithVal(lReplace));
    obj.locData.loc.(['ynmaligned_' obj.name])(idxWithVal(lReplace)) = newlocs.ynm(idxSubWithVal(lReplace));
    if isfield(newlocs,'znm')
        obj.locData.loc.(['znmaligned_' obj.name])(idxWithVal(lReplace)) = newlocs.znm(idxSubWithVal(lReplace));
    end
    obj.locData.loc.(['siteID_' obj.name])(idxWithVal(lReplace)) = obj.site.ID;
    obj.locData.loc.class(idxWithVal(lReplace)) = obj.site.ID;
    obj.locData.loc.(['rank_' obj.name])(idxWithVal(lReplace)) = obj.site.indList;
end
end