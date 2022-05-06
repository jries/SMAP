classdef LocMoFitGUI<interfaces.SEEvaluationProcessor
    properties
        fitter              % An SMLMModelFit object.
        numMod              % Number of component models.
        parsArgFieldnames   %
        lFnParsArgEdit      %
        fnParsArgColWidth   %
        layerFieldnames     %
        lFnLayerEdit        %
        currentLoadedModel  %
        alignSettings       % Converter for alignment.
        sourceModel         % The source function model of the image model.
    end
    methods
        function obj=LocMoFitGUI(varargin)
            obj@interfaces.SEEvaluationProcessor(varargin{:});
            addpath(genpath('./LocMoFit'))
            obj.propertiesToSave={'fitter', 'numMod', 'parsArgFieldnames', 'lFnParsArgEdit', 'fnParsArgColWidth', 'layerFieldnames', 'lFnLayerEdit', 'currentLoadedModel','sourceModel'};
            addlistener(obj, 'mParsArgModified', @mParsArgModified_callback);
        end
        
        
        function setGuiParameters(obj,p,setchildren,setmenu)
            obj.fitter = p.fitter;
            obj.fitter.updateVersion;
            initTabWhenLoading(obj);
            % For upgrade from SMLMModelFit to LocMoFit
            if ~isempty(p.anchorConvert.Data)
                col_source = p.anchorConvert.Data(:,1);
                col_source = replace(col_source, 'SMLMModelFitGUI', 'LocMoFitGUI');
                p.anchorConvert.Data(:,1) = col_source;
            end
            
            % Check the model type are consistent between the GUI and obj.
            m = 1;
            while isfield(p, ['modelType_' num2str(m)])
                ID = ['modelType_' num2str(m)];
                p.(ID).String = obj.guihandles.(ID).String;
                p.(ID).Value = obj.guihandles.(ID).Value;
                if iscell(p.(ID).String)
                    p.(ID).selection = p.(ID).String{p.(ID).Value};
                else
                    p.(ID).selection = p.(ID).String;
                end
                m = m+1;
            end
            setGuiParameters@interfaces.SEEvaluationProcessor(obj,p);
        end
        
        function initTabWhenLoading(obj)
            fitter = obj.fitter;
            numOfModel = fitter.numOfModel;
            inp = obj.getAllParameters;
            obj.setPar('loading',true);
            for m = 1:numOfModel
                modelnumberStr = num2str(m);
                if ~isfield(inp,['modelname_' modelnumberStr])
                    % hack the action of pressing '+' tab
                    eventData.EventName = 'SelectionChanged';
                    eventData.OldValue = [];
                    eventData.NewValue = [];
                    eventData.NewValue.Title = '+';
                    obj.guihandles.tabgroup.SelectionChangedFcn{1}(obj.guihandles.tabgroup,eventData,obj)
                end
                % load info to the GUI
                obj.currentLoadedModel = modelnumberStr;
                notify(obj.guihandles.(['modelload_' modelnumberStr]),'Action');
            end
            obj.setPar('loading',false);
        end
        
        function out=run(obj, inp, varargin)
            % varargin is used only when called by the run_addguitotab.mlx
            p = inputParser;
            addParameter(p,'onlySetUp',false, @islogical);
            addParameter(p,'forceDisplay',false, @islogical);
            addParameter(p,'keepParsVal',false, @islogical);
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
                fitter = obj.fitter;
                %% When pick site is checked, do nothing but pass on the old output
                for k = fitter.numOfModel:-1:1
                    lAllPicked(k) = inp.(['pickSite_' num2str(k)]);
                end
                lAnyPickedOn = any(lAllPicked);
                
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
                
                %% Deal with locs
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

                %% Alignment (for averaging)
                if inp.useAlignment
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
                %%
                %
                
                %% Display the result
                if inp.noFit||forceDisplay
                    if keepParsVal
                        tempAllParsArg = fitter.allParsArg; % this is to temporarily save the original allParsArg.
                    end
                    
                    out.parsInit = fitter.parsInit;
                    fitter.getDerivedPars;
                    
                    %% Get compensationFactor
                    % This factor compensate the number of locs between different channels
                    compensationFactor = zeros(size(locs.layer))';
                    for k = 1:fitter.numOfLayer
                        compensationFactor(locs.layer==k) = fitter.compensationFactor(k);
                    end
                    compensationFactor(~ismember(locs.layer, fitter.allModelLayer)) = [];
                    [~,~,parBestFit] = fitter.prepFit;
                    %     fitter.fitInfo.ELL = fitter.getELL(parBestFit', compensationFactor, 2);
                    
                    out.fitInfo = fitter.fitInfo;
                    out.externalInfo = fitter.externalInfo;
                end
                if ~lAnyPickedOn
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
            end
        end
        
        function makeGui(obj,varargin)
            %% init settings
            % create the SMLMModelFit obj
            % check the dim of the data and update the corresponding
            % setting in the fitter
            if isfield(obj.locData.loc,'znm')
                dataDim = 3;
            else
                dataDim = 2;
            end
            obj.fitter = SMLMModelFit('DataDim',dataDim);
            
            obj.currentLoadedModel = [];
            if ismac
                obj.guiPar.fontsize=12;
                obj.guiPar.Vrim=0;
                obj.guiPar.Vpos=3;
                obj.guiPar.FieldHeight=20;
            end
            
            %% ParArg table related things
            % field names in the parsArg table
            fnParsArg={'name','value','fix','lb','ub','type','min','max','label'};
            lFnParsArgEdit = logical([0 1 1 1 1 0 1 1 1]);
            fnParsArgColWidth = {50,35,20,30,30,30,30,30,30};
            obj.parsArgFieldnames = fnParsArg;
            obj.lFnParsArgEdit = lFnParsArgEdit;
            obj.fnParsArgColWidth = fnParsArgColWidth;
            % field names in the layer table
            fnLayer ={'layer','value','fix','lb','ub','min','max'};
            fnLayerColWidth = {30,35,20,30,30,30,30};
            lFnLayerEdit = logical([0 1 1 1 1 1 1]);
            obj.layerFieldnames = fnLayer;
            obj.lFnLayerEdit = lFnLayerEdit;
            
            if nargin >2 && varargin{2}==true %gui of tabs
                makeGui@interfaces.GuiModuleInterface(obj,varargin{1});
            else %real make GUI of main gui
                makeGui@interfaces.GuiModuleInterface(obj); %make the main GUI from the guidef definitions
                %Settings
                optimizernames=obj.fitter.supportedSolver;
                obj.guihandles.optimizer.String=optimizernames;
                % Use fminsearchbnd by default.
                obj.guihandles.optimizer.Value = find(strcmp(optimizernames, 'fminsearchbnd'));
                Data={};
                oldh=obj.guihandles.optimizerpar;
                pos=oldh.Position;
                htable=uitable(oldh.Parent,'Data',Data,'Position',pos);
                htable.ColumnEditable = true;
                htable.ColumnName = {'Parameters', 'Value'};
                htable.ColumnWidth = {100, 60};
                htable.CellSelectionCallback = {@optimizerSetting_selectionCallback,obj};
                htable.CellEditCallback = {@optimizerSetting_editCallback,obj};
                htable.RowName = [];
                delete(oldh);
                obj.guihandles.optimizerpar=htable;
                notify(obj.guihandles.optimizer, 'Action');
                
                %% Layer settings
                oldh = obj.guihandles.layerSetting;
                pos=oldh.Position;
                hLayer=uitable(oldh.Parent,'Data',{},'Position',pos);
                hLayer.ColumnName = fnLayer;
                hLayer.ColumnEditable = true;
                hLayer.CellEditCallback = {@layerSetting_callback,obj};
                hLayer.RowName = [];
                hLayer.ColumnWidth = fnLayerColWidth;
                delete(oldh);
                obj.guihandles.layerSetting=hLayer;
                
                %% add first module and + tab
                obj.guihandles.tabgroup=obj.guihandles.tab1.Parent;
                obj.addguitotab(1);
                obj.guihandles.addtab=uitab(obj.guihandles.tabgroup,'Title','+');
                obj.numMod = 1;                 % init of the model counts
                obj.guihandles.tabgroup.SelectionChangedFcn={@selectLayer_callback,obj};
                
                % Select the M1 tab by default since a user usually starts from loading a model.
                obj.guihandles.tabgroup.SelectedTab = obj.guihandles.tab1;
                
                %% Converter tab
                addconverttotab(obj);           % create the converter
                oldh=obj.guihandles.anchorConvert;
                
                % create the converter input table
                htable = createConvertTable(oldh, obj);
                obj.guihandles.anchorConvert=htable;
                
                
            end
            obj.setPar('initiated',true)
        end
        function pard=guidef(obj)
            if ismac
                dy=1;
            else
                dy=0;
            end
            pard.tab.tab1='Settings';
            
            pard.t_optimizer.object=struct('Style','text','String','Optimizer:');
            pard.t_optimizer.position=[1.5+dy,1];
            pard.t_optimizer.Width=1;
            pard.t_optimizer.tab='tab1';
            
            pard.optimizer.object=struct('Style','popupmenu','String',{{'text','t2'}}, 'Callback',{{@optimizer_callback,obj}});
            pard.optimizer.position=[1.5+dy,2];
            pard.optimizer.Width=1.5;
            pard.optimizer.tab='tab1';
            
            pard.loadfitting.object=struct('Style','pushbutton','String','load', 'Callback',{{@obj.load_callback}});
            pard.loadfitting.position=[1+dy,3.8];
            pard.loadfitting.Width=0.6;
            pard.loadfitting.Height=0.7;
            pard.loadfitting.tab='tab1';
            pard.loadfitting.TooltipString='Load previously saved settings.';
            
            pard.savefitting.object=struct('Style','pushbutton','String','save', 'Callback',{{@save_callback,obj}});
            pard.savefitting.position=[1+dy,4.4];
            pard.savefitting.Width=0.6;
            pard.savefitting.Height=0.7;
            pard.savefitting.tab='tab1';
            pard.savefitting.TooltipString='Save the current settings.';
            
            pard.setting_advanced.object=struct('Style','pushbutton','String','Advanced','Callback',{{@setting_advanced_callback,obj}});
            pard.setting_advanced.position=[1.8+dy,3.8];
            pard.setting_advanced.Width=1.2;
            pard.setting_advanced.Height=0.7;
            pard.setting_advanced.tab='tab1';
            pard.setting_advanced.TooltipString='Advanced settings.';
            
            pard.noFit.object=struct('Style','checkbox','String','Review only','Value',0);
            pard.noFit.position=[3.5+dy,3.5];
            pard.noFit.Width=1.3;
            pard.noFit.Height=1;
            pard.noFit.tab='tab1';
                       
            pard.useAlignment.object=struct('Style','checkbox','String','Transform','Value',0);
            pard.useAlignment.position=[4.5+dy,3.5];
            pard.useAlignment.Width=1.5;
            pard.useAlignment.Height=1;
            pard.useAlignment.tab='tab1';
            pard.useAlignment.TooltipString = 'Perform the transformation of the site based on the model.';
            
            pard.setting_alignment.object=struct('Style','pushbutton','String','...','Callback',{{@setting_alignment_callback,obj}});
            pard.setting_alignment.position=[4.3+dy,4.7];
            pard.setting_alignment.Width=0.25;
            pard.setting_alignment.Height=0.5;
            pard.setting_alignment.tab='tab1';
                                  
            pard.optimizerpar.object=struct('Style','text','String','');
            pard.optimizerpar.position=[4.5+dy,1];
            pard.optimizerpar.Width=2.5;
            pard.optimizerpar.Height=3;
            pard.optimizerpar.tab='tab1';
            
            pard.addRowOptimizer.object=struct('Style','pushbutton','String','+', 'Callback',{{@addRowOptimizer_callback,obj}});
            pard.addRowOptimizer.position=[5+dy,3];
            pard.addRowOptimizer.Width=0.2;
            pard.addRowOptimizer.Height=0.4;
            pard.addRowOptimizer.tab='tab1';
            
            pard.rmRowOptimizer.object=struct('Style','pushbutton','String','-', 'Callback',{{@rmRowOptimizer_callback,obj}});
            pard.rmRowOptimizer.position=[5+dy,3.2];
            pard.rmRowOptimizer.Width=0.2;
            pard.rmRowOptimizer.Height=0.4;
            pard.rmRowOptimizer.tab='tab1';
            
            pard.t_layerSetting.object=struct('Style','text','String','Layer background:');
            pard.t_layerSetting.position=[6+dy,1];
            pard.t_layerSetting.Width=1.5;
            pard.t_layerSetting.tab='tab1';
            
            pard.layerSetting.object=struct('Style','text','String','');
            pard.layerSetting.position=[10+dy,1];
            pard.layerSetting.Width=3.9;
            pard.layerSetting.Height=4;
            pard.layerSetting.tab='tab1';
        end
        
        function addguitotab(obj,number)
            %% RUN_ADDGUITOTAB Adding a model (number) to the tab group.
            %
            tag=['M' num2str(number)];
            obj.guihandles.(['tab' num2str(number)])=uitab(obj.guihandles.tabgroup,'Title',tag,'Tag',tag);
            guidefhere=addnumbertofield(guidefmodel(obj,number),number);
            Vrimold=obj.guiPar.Vrim;handleold=obj.handle;
            obj.guiPar.Vrim=0;
            if ismac
                obj.guiPar.FieldHeight=20;
                obj.guiPar.fontsize=12;
                obj.guiPar.Xrim=0;
                tabh=[50 160];
                
            else
                tabh=[20 120];
            end
            obj.handle=obj.guihandles.(['tab' num2str(number)]);
            fitter = obj.fitter;
            obj.makeGui(guidefhere,1);
            obj.fitter = fitter;
            obj.fitter.linkedGUI = obj;
            obj.handle=handleold;
            obj.guiPar.Vrim=Vrimold;
            %initialize parameter table. Here change the looks!
            hpar=obj.guihandles.(['tabpar_' num2str(number)]);
            hsettings=obj.guihandles.(['tabsettings_' num2str(number)]);
            
            % parameter tab
            
            % flatten the guihandles
            obj.guihandles.(['partable_' num2str(number)])=uitable(hpar,'Data',{});
            obj.guihandles.(['partable_' num2str(number)]).ColumnName = obj.parsArgFieldnames;
            obj.guihandles.(['partable_' num2str(number)]).Position(3:4)=obj.guihandles.(['partable_' num2str(number)]).Position(3:4)-tabh;
            obj.guihandles.(['partable_' num2str(number)]).ColumnEditable = obj.lFnParsArgEdit;
            obj.guihandles.(['partable_' num2str(number)]).RowName = [];
            
            % button for saving the model
            obj.guihandles.(['savePar_' num2str(number)])=uicontrol(hpar,'Style','pushbutton','String','Save');
            obj.guihandles.(['savePar_' num2str(number)]).Position = [20 0 40 20];
            obj.guihandles.(['savePar_' num2str(number)]).Callback = {@savePar_callback, obj, number};
            
            % button for loading the model
            obj.guihandles.(['loadPar_' num2str(number)])=uicontrol(hpar,'Style','pushbutton','String','Load');
            obj.guihandles.(['loadPar_' num2str(number)]).Position = [60 0 40 20];
            obj.guihandles.(['loadPar_' num2str(number)]).Callback = {@loadPar_callback, obj, number};
            
            % toggle button for selecting a site
            obj.guihandles.(['pickSite_' num2str(number)])=uicontrol(hpar,'Style','togglebutton','String','Pick site');
            obj.guihandles.(['pickSite_' num2str(number)]).Position = [180 0 60 20];
            
            % button for previewing the model
            obj.guihandles.(['forceDisplay_' num2str(number)])=uicontrol(hpar,'Style','pushbutton','String','Preview');
            obj.guihandles.(['forceDisplay_' num2str(number)]).Position = [240 0 60 20];
            obj.guihandles.(['forceDisplay_' num2str(number)]).Callback = {@displayModel_callback, obj};
            
            % internal settings tab
            obj.guihandles.(['settingstable_' num2str(number)])=uitable(hsettings,'Data',{});
            obj.guihandles.(['settingstable_' num2str(number)]).ColumnName = {'setting','value'};
            obj.guihandles.(['settingstable_' num2str(number)]).Position(3:4)=obj.guihandles.(['settingstable_' num2str(number)]).Position(3:4)-tabh;
            obj.guihandles.(['settingstable_' num2str(number)]).ColumnEditable = true;
            obj.guihandles.(['settingstable_' num2str(number)]).RowName = [];
            obj.guihandles.(['settingstable_' num2str(number)]).CellSelectionCallback = {@settingTable_cellSelectionCallback, obj, num2str(number)};
            
            % Buttons for adding removing a row of setting
            obj.guihandles.(['addInternalSetting_' num2str(number)])=uicontrol(hsettings,'Style','pushbutton','String','+');
            obj.guihandles.(['addInternalSetting_' num2str(number)]).Position = [260 0 20 20];
            obj.guihandles.(['addInternalSetting_' num2str(number)]).Callback = {@addInternalSetting_callback, obj, num2str(number)};
            obj.guihandles.(['rmInternalSetting_' num2str(number)])=uicontrol(hsettings,'Style','pushbutton','String','-');
            obj.guihandles.(['rmInternalSetting_' num2str(number)]).Position = [280 0 20 20];
            obj.guihandles.(['rmInternalSetting_' num2str(number)]).Callback = {@rmInternalSetting_callback, obj, num2str(number)};
        end
        
        function addconverttotab(obj)
            tag='Convert';
            obj.guihandles.converter=uitab(obj.guihandles.tabgroup,'Title',tag,'Tag',tag);
            guidefhere=guidefconvert(obj);
            Vrimold=obj.guiPar.Vrim;handleold=obj.handle;
            obj.guiPar.Vrim=0;
            obj.handle=obj.guihandles.converter;
            obj.makeGui(guidefhere,1);
            obj.handle=handleold;
            obj.guiPar.Vrim=Vrimold;
        end
        
        function parId = loadParTable(obj, htable, fitter, modelnumber)
            % get parId and update the GUIParTable
            [parId,subParsArgTemp] = fitter.getAllParId(modelnumber, 'form', 'long');
            htable.Data = struct2Data(subParsArgTemp);
            htable.CellEditCallback = {@parSetting_callback,obj, modelnumber};
            htable.ColumnEditable = obj.lFnParsArgEdit;
            htable.ColumnWidth = obj.fnParsArgColWidth;
        end
        
        %%
        function mParsArgModified_callback(obj,b)
            modelbnumberStr = num2str(b.modelID);
            htable = obj.guihandles.(['partable_' modelbnumberStr ]);
            parId = obj.loadParTable(htable, obj.fitter, b.modelID);
            % here the target options in the converter table are defined.
            hConvert = obj.guihandles.anchorConvert;
            optionTarget = unique([hConvert.ColumnFormat{3} parId]);
            hConvert.ColumnFormat{3} = optionTarget;
            obj.guihandles.anchorConvert=hConvert;
        end
        
        %% advanced settings
        function setAdvanceSetting(obj, par, value)
            obj.advanceSetting.(par) = value;
        end
        
        function set.fitter(obj,value)
            obj.fitter = value;
            if isfield(obj.P.par.mainGui.content.children.guiSites.children.Segment.processors,'SimulateSites')
                simulateSites = obj.P.par.mainGui.content.children.guiSites.children.Segment.processors.SimulateSites;
                if strcmp(simulateSites.guihandles.useFitter_button.Visible, 'off')
                    obj.P.par.mainGui.content.children.guiSites.children.Segment.processors.SimulateSites.initGui;
                end
            end
        end
        %% Updating the layer settings for the GUI
        function updateLayer(obj)
            % Update the layers in use.
            % Get the fitter.
            fitter = obj.fitter;
            hLayer=obj.guihandles.layerSetting; % Hande of the layer table.
            % Get the current settings for the layers.
            layerParsArg = fitter.subParsArg([]);
            % Convert cellstr into chr array so that uitalble can handle.
            fn = obj.layerFieldnames;
            layerStr = num2str(layerParsArg.model);
            layerParsArgDisp.layer = char(strcat('layer', layerStr(:,2)));
            for k = 2:length(fn)
                if iscell(layerParsArg.(fn{k}))
                    layerParsArgDisp.(fn{k}) = char(layerParsArg.(fn{k}));
                else
                    layerParsArgDisp.(fn{k}) = layerParsArg.(fn{k});
                end
            end
            hLayer.Data = struct2Data(layerParsArgDisp);
            hLayer.ColumnEditable = obj.lFnLayerEdit;
            % save the table back to the GUI
            obj.guihandles.layerSetting = hLayer;
            layerParId = {};
            for k = length(layerParsArg.model):-1:1
                layerParId(k) = fitter.getAllParId(layerParsArg.model(k));
            end
            %% update the convert table
            hConvert = obj.guihandles.anchorConvert;
            colFormat = hConvert.ColumnFormat{3};
            % remove old layer parID
            lLayerParID = contains(colFormat,'.offset.weight');
            if any(lLayerParID)
                colFormat(lLayerParID) = [];
            end
            % add new layer parID
            optionTarget = unique([colFormat layerParId{:}]);
            hConvert.ColumnFormat{3} = optionTarget;
            obj.guihandles.anchorConvert=hConvert;
        end
        % callback for loading
        function load_callback(obj,a,b,path2setting)
            % get the path of the saved info
            if obj.getPar('initiated')
                if ~exist('path2setting', 'var')
                    [file,path2setting] = uigetfile({'*_SMLMModelFit.mat;*_LocMoFit.mat','LocMoFit settings'}, 'Select an SMLMModelFit file...', '');
                else
                    [path2setting,file,ext] = fileparts(path2setting);
                    file = [file ext];
                    path2setting = [path2setting filesep];
                end
                temp = load([path2setting file]);
                obj.setGuiParameters(temp.par);
                obj.fitter = temp.par.fitter;
                obj.fitter.linkedGUI = obj;
                obj.fitter.loadListener;
                for k = 1:obj.fitter.numOfModel
                    obj.fitter.model{k}.loadListener;
                end
            end
        end

        function updateGUI_fromLocMoFitObj(obj)
            fitter = obj.fitter;
            for m = 1:fitter.numOfModel
                htable = obj.guihandles.(['partable_' num2str(m)]);
                obj.loadParTable(htable, fitter, m);
                obj.updateAdvanceTab(fitter, m);
            end
        end

        function updateGUI_convert_fromLocMoFitObj(obj)
            % Uupdate GUI: update the convert table based on the LocMoFit
            % Obj.
            % Known issue: targets start with 'usr_' cannot be processed.
            fitter = obj.fitter;
            allSourceFitter = fitter.converterSource;
            loaded_beforeMe = obj.find_LocMoFit_beforeMe;
            fn = fieldnames(loaded_beforeMe);
            fitterGUI_order = {};
            for k = 2:length(allSourceFitter)
                for m = 1:length(loaded_beforeMe)
                    l = isequal(allSourceFitter{k}, loaded_beforeMe.(fn{m}));
                    if l
                        fitterGUI_order{k} = fn{m};
                    end
                end
            end
            numOfRules = length(fitter.converterRules.target);
            data = cell(numOfRules,4);
            data(:,1:3) = [fitterGUI_order(fitter.converterRules.target_Id)' fitter.converterRules.rule_raw' fitter.converterRules.target'];
            obj.guihandles.anchorConvert.Data = data;
        end

        function updateAdvanceTab(obj,fitter, modelnumberStr)
            m = modelnumberStr;
            htable = obj.guihandles.(['settingstable_' num2str(m)]);
            internalPar_list = fitter.getModelInternalSettingList(m);
            data = {};
            for k = length(internalPar_list):-1:1
                currentPar = internalPar_list{k};
                val = fitter.getModelInternalSetting(m, currentPar);
                data(k,:) = {currentPar val};
            end
            htable.Data = data;
        end

        function fitter_stack = find_LocMoFit_beforeMe(obj)
            se = obj.locData.SE;
            eval_ = se.processors.eval;
            loadedPlugins = fieldnames(eval_.children);
            loadedLocMoFitGUI = loadedPlugins(startsWith(loadedPlugins, 'LocMoFit'));
            currentName = obj.name;
            ind_current = find(ismember(loadedLocMoFitGUI, currentName));
            if ind_current > 1
                loadedLocMoFitGUI_beforeMe = loadedLocMoFitGUI(1:ind_current-1);
                fn = fieldnames(se.processors.eval.children);
                for k = 1:length(loadedLocMoFitGUI_beforeMe)
                    fitter_stack.(loadedLocMoFitGUI_beforeMe{k}) = se.processors.eval.children.(loadedLocMoFitGUI_beforeMe{k}).fitter;
                end
            else
                fitter_stack = [];
            end
            
        end
    end
    events
        mParsArgModified
    end
end

%% Callbacks
function parSetting_callback(a,b,obj, modelnumber)
% Callback for editing parsArg of the model
% First get the fitter obj, and then reset the initial guess based on what
% are saved. This is to prevent any interference between different sites.
% Next detect the position of the modified value, check its counterpart in
% fitter.allParsArg, and overwrite it. Finally the initial guess is saved
% and the fitter of the obj replaced.

fitter = obj.fitter;
fitter.resetInit;
fn = obj.parsArgFieldnames;
indices = b.Indices;
data = a.Data;
name = strtrim(data{indices(1),strcmp(fn,'name')});
type = strtrim(data{indices(1),strcmp(fn,'type')});
[~,idx] = fitter.wherePar(['pars.m' num2str(modelnumber) '.' type '.' name]);
fitter.allParsArg.(fn{indices(2)})(idx) = b.NewData;
fitter.saveInit;
obj.fitter = fitter;
end

function selectLayer_callback(tabgroup,eventdata,obj)
% if + tab selected this makes a new model
layertitle=(eventdata.NewValue.Title);
if strcmp(layertitle,'+')
    numtabs=length(tabgroup.Children);
    obj.addguitotab(numtabs-2)
    obj.numMod = numtabs-2;   % save the number of model
    s=1:length(tabgroup.Children);
    % shift the order of table
    s(end-2)=s(end);
    s(end)=s(end)-2;
    s(end-1)=s(end);
    s(end)=s(end)+1;
    tabgroup.Children=tabgroup.Children(s);
    tabgroup.SelectedTab=tabgroup.Children(end-2);
end
end

function layerSetting_callback(a,b,obj)
fitter = obj.fitter;
fn=obj.layerFieldnames;
indices = b.Indices;
layer = a.Data{indices(1),1};
layer = strrep(layer,'layer','9');
offsetForm = obj.fitter.getAdvanceSetting(['m' num2str(layer) '_background']);
[~,idx] = fitter.wherePar(['pars.m' layer '.offset.' offsetForm]);
fitter.allParsArg.(fn{indices(2)})(idx) = b.NewData;
obj.fitter = fitter;
end

function optimizerSetting_selectionCallback(a,b,obj)
try
    selectedRow = b.Indices(1);
    obj.setPar(['selectedRowOptimizer_' obj.name], selectedRow)
catch
end
end

function optimizerSetting_editCallback(a,b,obj)
obj.setPar(['lOptimizerSettingEdited_' obj.name], true)
end

%% Format conversion
% From uitableData to parsArg struct
function parsArg = data2ParsArg(uitableData,fn,model)
parsArgOri = cell2struct(uitableData,fn(2:end),2);
for k = 2:length(fn)
    if ischar([parsArgOri.(fn{k})])
        tempField = strtrim({parsArgOri.(fn{k})});
    else
        tempField = [parsArgOri.(fn{k})];
    end
    parsArg.(fn{k}) = tempField';
end
parsArg.model = repelem(model, length(parsArg.(fn{2})))';
parsArg = orderfields(parsArg, fn);
end

%% guidef related


% Call back for saving SMLMModelFit settings

function save_callback(a,b,obj)
[file,path] = uiputfile('*_LocMoFit.mat', 'Save as', '');
obj.fitter.resetInit;
gm = obj.getPar('mainGui');
ge = gm.children.guiSites.children.eval;
pe = ge.getGuiParameters(true);
par = pe.children.(obj.name);
var2save = {'par'};
save([path file], [var2save{:}],'-v7.3');
end

%         Create a new row in the optimizer table

function addRowOptimizer_callback(a,b,obj)
htable = obj.guihandles.optimizerpar;
htable.Data = [htable.Data; {[],[]}];
obj.guihandles.optimizerpar=htable;
end

%         Delete a row in the optimizer table

function rmRowOptimizer_callback(a,b,obj)
htable = obj.guihandles.optimizerpar;
data = htable.Data;
selectedRow = obj.getPar(['selectedRowOptimizer_'  obj.name]);
fitter = obj.fitter;
try
    data(selectedRow,:) = [];
catch
end
htable.Data = data;
obj.guihandles.optimizerpar=htable;
end

%         Call back for selecting the optimizer

function optimizer_callback(a,b,obj)
% !!! later here should update the options for the optimizer settings
% based on the optimizer selected.
selectedOptimizer = b.Source.String{b.Source.Value};
hSetting = obj.guihandles.optimizerpar;
hSetting.ColumnFormat = {getSolverOption(selectedOptimizer),[]};
obj.setPar(['lOptimizerSettingEdited_' obj.name], true);
end

function setting_advanced_callback(a,b,obj)
fig = figure('Name','Advance setting');
advanceSetting = obj.fitter.advanceSetting;
fn = fieldnames(advanceSetting);
nrow = length(fn);
fig.Position(3:4) = [320 nrow*22+60];
guihandles = [];
for k = 1:nrow
    option = advanceSetting.(fn{k}).option;
    value = advanceSetting.(fn{k}).value;
    name = advanceSetting.(fn{k}).name;
    t = ['uit_' num2str(k)];
    ui = ['ui_' num2str(k)];
    lineHeight = 0.73;
    if islogical(value)
        guihandles.(t) = uicontrol(fig,'Style','checkbox','String',name,'Value',value);
        guihandles.(t).Position = [1.2 (k-1)*lineHeight 1.5 0.9];
    else
        if length(option)>1
            idxVal = find(strcmp(value,option));
            guihandles.(t) = uicontrol(fig,'Style','text','String',name);
            guihandles.(t).Position = [1.2 (k-1)*lineHeight 1.5 0.9];
            guihandles.(ui) = uicontrol(fig,'Style','popupmenu','Value',idxVal,'String',option);
            guihandles.(ui).Position = [2.7 (k-1)*lineHeight 1 0.9];
        else
            guihandles.(t) = uicontrol(fig,'Style','text','String',name);
            guihandles.(t).Position = [1.2 (k-1)*lineHeight 1.5 0.9];
            guihandles.(ui) = uicontrol(fig,'Style','edit','String',num2str(value));
            guihandles.(ui).Max = 2;
            guihandles.(ui).Position = [2.7 (k-1)*lineHeight 1 0.9];
        end
    end
    
    if islogical(value)
        obj.guihandles.([tag 't_' fn{k}]) = guihandles.(t);
    end
    
    obj.guihandles.(fn{k}) = guihandles.(ui);
end

guihandles.save = uicontrol(fig,'Style','pushbutton','String','Save','Callback',{@saveAdvanceSettings_callback, fig, obj});
guihandles.save.Position = [3.2 k*lineHeight+0.3 0.5 0.9];

guiStyle(guihandles, fieldnames(guihandles),'FieldHeight',20);
end

function saveAdvanceSettings_callback(a,b,fig,obj)
fitter = obj.fitter;
advanceSetting = fitter.advanceSetting;
fn = fieldnames(advanceSetting);
for k = 1:length(fn)
    oneCtrl = obj.guihandles.(fn{k});
    switch oneCtrl.Style
        case 'popupmenu'
            fitter.setAdvanceSetting(fn{k},oneCtrl.String{oneCtrl.Value});
        case 'edit'
            fitter.setAdvanceSetting(fn{k},str2num(oneCtrl.String));
        case 'checkbox'
            fitter.setAdvanceSetting(fn{k},oneCtrl.String);
    end
end
try
    obj.updateLayer
catch
    warning('Layer(s) is not updated.')
end
close(fig)
end

%         Call back for alignment settings

function setting_alignment_callback(a,b,obj)
fig = figure(514);
fig.Name = 'Transformation settings';
fig.Position(3:4) = [400 360];
guihandles.uit = uitable(fig);
guihandles.uit = createConvertTable(guihandles.uit, fig);
guihandles.uit.CellEditCallback = {@alignSettingEdit_callback,3,obj};

% Update the convert tab.
parId = {'post_x', 'post_y', 'post_z', 'post_scale', 'post_zrot'};
optionTarget = unique([guihandles.uit.ColumnFormat{2} parId]);
guihandles.uit.ColumnFormat{2} = optionTarget;
guihandles.uit.Data = obj.alignSettings;
guihandles.uit.ColumnWidth = {100 100 100};



guihandles.button_add = uicontrol('Style','pushbutton','String','+','Callback',{@addNewRuleAlign_callback,fig,obj});
guihandles.button_add.Position = [1 10 0.3 0.8];

guihandles.button_rm = uicontrol('Style','pushbutton','String','-','Callback',{@rmRuleAlign_callback,fig,obj});
guihandles.button_rm.Position = [1.3 10 0.3 0.8];
guihandles.uit.Position = [1 0 3.7 10];
guiStyle(guihandles, fieldnames(guihandles));
end

function addNewRuleAlign_callback(a,b,fig,obj)
htable = findobj(fig, 'Type', 'uitable');
htable.Data = [htable.Data; {[],[],[]}];
end

function rmRuleAlign_callback(a,b,fig,obj)
htable = findobj(fig, 'Type', 'uitable');
data = htable.Data;
temp = findobj(fig,'Type','uicontrol','-and','String','selectedRowConvert');
try
    selectedRow = temp.Value;
    data(selectedRow,:) = [];
    htable.Data = data;
    obj.alignSettings = htable.Data;
catch
    display('Please select a row first.')
end

end

function alignSettingEdit_callback(a,b,k,obj)
obj.alignSettings = a.Data;
end

%% addguitotab related
function out=addnumbertofield(in,number)
% A helper function. obj.guihandles needs to be flat. Add number to field name.
fn=fieldnames(in);
for k=1:length(fn)
    if strcmp(fn{k},'tab')
        out.tab=addnumbertofield(in.(fn{k}),number);
    else
        out.([fn{k} '_' num2str(number)])=in.(fn{k});
    end
end
end

% GUI for the model tab
function pard=guidefmodel(obj,number)
if ismac
    dy=1.5;
else
    dy=0;
end
pard.tab.tabmodel='Model';
pard.modelname.object=struct('Style','edit','String','');
pard.modelname.position=[3+dy,1];
pard.modelname.Width=2;
pard.modelname.tab=['tabmodel_' num2str(number)];

pard.modelload.object=struct('Style','pushbutton','String','load model','Callback',{{@loadmodel_callback,obj}});
pard.modelload.position=[3+dy,3];
pard.modelload.Width=1;
pard.modelload.tab=['tabmodel_' num2str(number)];

pard.modelType.object=struct('Style','popupmenu','String',{'Type'},'value', 1,'Callback',{{@modType_callback,obj,number}});
pard.modelType.position=[3+dy,4];
pard.modelType.Width=1;
pard.modelType.tab=['tabmodel_' num2str(number)];
pard.modelType.Tooltip='Convert the model to the type you specify.';
pard.modelType.Enable = 'off';

pard.layert.object=struct('Style','text','String','Layer');
pard.layert.position=[4+dy,1];
pard.layert.Width=2;
pard.layert.tab=['tabmodel_' num2str(number)];
pard.layert.Tooltip = 'The layer that the model is fit to.';

pard.layer.object=struct('Style','popupmenu','String',{{'1','2','3','4','5','6'}}, 'Callback', {{@layer_callback,obj,number}});
pard.layer.position=[4+dy,3.5];
pard.layer.Width=1;
pard.layer.tab=['tabmodel_' num2str(number)];
pard.layer.Tooltip = "This corresponds to the layer defined in the 'Render' tab.";
pard.layer.Enable = 'off';

pard.pixelsizefitt.object=struct('Style','text','String','Pixel size');
pard.pixelsizefitt.position=[5+dy,1];
pard.pixelsizefitt.Width=2;
pard.pixelsizefitt.tab=['tabmodel_' num2str(number)];
pard.pixelsizefitt.Tooltip = 'The bin size for binning the model.';

pard.pixelsizefit.object=struct('Style','edit','String','10','Callback',{{@modelSettingsEdited_callback,obj}});
pard.pixelsizefit.position=[5+dy,3.5];
pard.pixelsizefit.Width=1;
pard.pixelsizefit.tab=['tabmodel_' num2str(number)];
pard.pixelsizefit.Tooltip = 'The smaller the more precise.';
pard.pixelsizefit.Enable = 'off';

pard.sigma_fitt.object=struct('Style','text','String','Sigma of Gaussian filtering');
pard.sigma_fitt.position=[6+dy,1];
pard.sigma_fitt.Width=2;
pard.sigma_fitt.tab=['tabmodel_' num2str(number)];
pard.sigma_fitt.Tooltip = 'The Gaussian filtering makes the model smoother.';

pard.sigma_fit.object=struct('Style','edit','String','15','Callback',{{@modelSettingsEdited_callback,obj}});
pard.sigma_fit.position=[6+dy,3.5];
pard.sigma_fit.Width=1;
pard.sigma_fit.tab=['tabmodel_' num2str(number)];
pard.sigma_fit.Enable = 'off';

pard.sigmaOrFactor.object=struct('Style','text','String','');
pard.sigmaOrFactor.position=[6+dy,4.5];
pard.sigmaOrFactor.Width=1;
pard.sigmaOrFactor.Height=2;
pard.sigmaOrFactor.tab=['tabmodel_' num2str(number)];
pard.sigmaOrFactor.Tooltip = 'Select the mode of sigma to use.';

pard.useSigma.object=struct('Style','radiobutton', 'Value', 1,'Callback',{{@modelSettingsEdited_callback,obj}});
pard.useSigma.position=[1+dy,1];
pard.useSigma.Width=1;
pard.useSigma.tab=['tabmodel_' num2str(number)];
pard.useSigma.Visible = 'off';

pard.sigmaFactor_fitt.object=struct('Style','text','String','Factor of Gaussian filtering');
pard.sigmaFactor_fitt.position=[7+dy,1];
pard.sigmaFactor_fitt.Width=2;
pard.sigmaFactor_fitt.tab=['tabmodel_' num2str(number)];
pard.sigmaFactor_fitt.Tooltip = 'The Gaussian filtering makes the model smoother.';

pard.useSigmaFactor.object=struct('Style','radiobutton','Value', 0,'Callback',{{@modelSettingsEdited_callback,obj}});
pard.useSigmaFactor.position=[2+dy,1];
pard.useSigmaFactor.Width=1;
pard.useSigmaFactor.tab=['tabmodel_' num2str(number)];
pard.useSigmaFactor.Visible = 'off';

pard.sigmaFactor_fit.object=struct('Style','edit','String','1','Callback',{{@modelSettingsEdited_callback,obj}});
pard.sigmaFactor_fit.position=[7+dy,3.5];
pard.sigmaFactor_fit.Width=1;
pard.sigmaFactor_fit.tab=['tabmodel_' num2str(number)];
pard.sigmaFactor_fit.Enable = 'off';

% the sub tabs for the model are created here
pard.tab.tabpar='Parameters';
pard.tab.tabsettings='Advance';
end

% Actions when loading a model
function loadmodel_callback(a,b,obj)
% Executed when loading a model.
% Load model based on the model number (ID).
currentLoadedModel = obj.currentLoadedModel;
if isempty(currentLoadedModel)
    modelnumber=(a.Parent.Parent.Parent.Title(2:end)); %hack to get the right tab
    filter = {'*.m;*.mlx;*_img.mat;*.png','Supported formats (**.m,*.mlx,*_img.mat,*.png)';...
        '*.m;*.mlx','Code files (*.m,*.mlx)'; ...
        '*_img.mat','Matrix files (*_img.mat)'; ...
        '*.png','Image files (*.png)'; ...
        };
    fnold=obj.guihandles.(['modelname_' modelnumber]).String;
    if isempty(fnold)
        [~,fnold] = LocMoFit.getModelList;
    end
    [f,p]=uigetfile(filter,'Select a geometric model',fnold);
    if ~f %no model selected: return
        disp('Cancelled by the user. No model loaded.')
        return
    end
    obj.guihandles.(['modelname_' modelnumber]).String=[p f];
else
    modelnumber = currentLoadedModel;
end
% If the model is loaded from a saved file, then skip addModel and load the
% model obj from the file, otherwise create a new obj.
if obj.getPar('loading')
    initmodel(obj, modelnumber,'skipAddModel',true);
else
    initmodel(obj, modelnumber);
end

% Update the model type options.
switch obj.fitter.model{str2num(modelnumber)}.modelType
    case 'image'
        if isempty(obj.sourceModel)
            modelOption ={'image'};
        else
            modelOption ={obj.sourceModel{str2num(modelnumber)}.modelObj.modelTypeOption{:} 'image'};
        end
    case 'locsImg'
        modelOption ={obj.fitter.model{str2num(modelnumber)}.modelObj.modelTypeOption{:}};
    otherwise
        modelOption ={obj.fitter.model{str2num(modelnumber)}.modelObj.modelTypeOption{:} 'image'};
end
indType = find(strcmp(modelOption,obj.fitter.model{str2num(modelnumber)}.modelType));
obj.guihandles.(['modelType_' modelnumber]).String = modelOption;
obj.guihandles.(['modelType_' modelnumber]).Value = indType;

% Update the layers in use.
obj.updateLayer;
end

%    Initiate the model
function initmodel(obj, modelnumber,varargin)
% Initiate the model.
% Initiate the function.
p = inputParser;
addParameter(p, 'skipAddModel',false)
parse(p,varargin{:});
results = p.Results;
fitter = obj.fitter;
% Get model's path.
modPath = obj.guihandles.(['modelname_' modelnumber]).String;
if ~strcmp(modPath,'Loaded')&&~results.skipAddModel
    
    %% Decide which subclass of the SMLMModel to use based on the input model.
    [filePath, modelFun, ext] = fileparts(obj.guihandles.(['modelname_' modelnumber]).String);
    
    % Add the root of the model script as a path
    
    %     modelFun = str2func(modelFun);
    %     tempModelObj = modelFun();
    if ismember(ext,{'.png', '.bmp', '.tif', '.tiff'}) || (strcmp(ext,{'.mat'}) && endsWith(modelFun,'_img'))
        if strcmp(ext, '.mat')
            img = load(modPath,'img');
            img = img.img;
        else
            img = imread(modPath);
        end
        geoModeltemp = imageModel(img);
    else
        addpath(filePath);
        geoModeltemp = functionModel(modPath);
    end
    
    if isempty(obj.locData.loc)
        % When there is not localizations, assuming the user is loading the
        % geometric model for simulations.
        fitter.dataDim = geoModeltemp.dimension;
    end

    % If the position of the model has been occupied, replace the occupying model.
    if ~isempty(fitter.model)&&~length(fitter.model)<str2double(modelnumber)
        fitter.changeModel(geoModeltemp, str2double(modelnumber));
    else
        fitter.addModel(geoModeltemp);
    end
end
% Set up the GUI of the model tab.
enableSettings(obj,modelnumber,fitter);
modelnumberStr = modelnumber;
modelnumber = str2double(modelnumber);
% Load model's info to the GUI.
obj.guihandles.(['layer_' modelnumberStr]).Value = fitter.model{modelnumber}.layer;
if ismember(fitter.model{modelnumber}.modelType,{'discrete','continuous'})
    obj.guihandles.(['sigma_fit_' modelnumberStr]).String = fitter.model{modelnumber}.sigma;
    obj.guihandles.(['sigmaFactor_fit_' modelnumberStr]).String = num2str(fitter.model{modelnumber}.sigmaFactor);
end
obj.guihandles.(['pixelsizefit_' modelnumberStr]).String = fitter.model{modelnumber}.pixelSize;

%             _for the parameter table_
% Update the model parameter table when the model has been loaded.
htable=obj.guihandles.(['partable_' modelnumberStr]); %handle of table
parId = obj.loadParTable(htable, fitter, modelnumber);
obj.fitter = fitter;

%             _for the internal setting table (in the advanced tab)_
% Assign the callback function and default options.
hSettingTable=obj.guihandles.(['settingstable_' modelnumberStr]); %handle of table
hSettingTable.CellEditCallback = {@settingTable_cellEditCallback,obj, modelnumber};

if isempty(fitter.getModelInternalSettingList(modelnumber))
    hSettingTable.ColumnFormat{1} = [];
else
    % !!!200528: I was asked to input a row vector instead of column. Keep this in
    % mind if things happend again.
    hSettingTable.ColumnFormat{1} = cellstr(fitter.getModelInternalSettingList(modelnumber))';
end

% Update the convert tab.
hConvert = obj.guihandles.anchorConvert;
optionTarget = unique([hConvert.ColumnFormat{3} parId]);
hConvert.ColumnFormat{3} = optionTarget;
obj.guihandles.anchorConvert=hConvert;
end

% Enable settings
function enableSettings(obj,modelnumber,fitter)
% Enable settings according to the model type.
modelType = fitter.model{str2double(modelnumber)}.modelType;
if isequal(modelType, 'image')
    % for image models
    obj.guihandles.(['layer_' modelnumber]).Enable = 'on';
    obj.guihandles.(['pixelsizefit_' modelnumber]).Enable = 'on';
    obj.guihandles.(['sigma_fit_' modelnumber]).Enable = 'off';
    obj.guihandles.(['sigmaFactor_fit_' modelnumber]).Enable = 'off';
elseif isequal(modelType, 'continuous')||isequal(modelType, 'background')
    % for continuous models
    obj.guihandles.(['layer_' modelnumber]).Enable = 'on';
    obj.guihandles.(['pixelsizefit_' modelnumber]).Enable = 'on';
    obj.guihandles.(['sigma_fit_' modelnumber]).Enable = 'on';
    obj.guihandles.(['sigmaFactor_fit_' modelnumber]).Enable = 'off';
else
    % for discrete/discretized models
    obj.guihandles.(['layer_' modelnumber]).Enable = 'on';
    obj.guihandles.(['pixelsizefit_' modelnumber]).Enable = 'off';
    obj.guihandles.(['pixelsizefit_' modelnumber]).String = 1;
    obj.guihandles.(['sigma_fit_' modelnumber]).Enable = 'off';
    obj.guihandles.(['sigmaFactor_fit_' modelnumber]).Enable = 'on';
    
    % for point models, user can choose to set sigma or sigma factor
    % here create the UI for the selection based on radiobutton
    obj.guihandles.(['useSigmaFactor_' modelnumber]).Visible = 'on';
    obj.guihandles.(['useSigma_' modelnumber]).Visible = 'on';
    hParent = obj.guihandles.(['sigmaOrFactor_' modelnumber]).Parent;
    obj.guihandles.(['sigmaOrFactor_' modelnumber]) = uibuttongroup('BorderType','none');
    obj.guihandles.(['sigmaOrFactor_' modelnumber]).SelectionChangedFcn = {@sigmaOrFactor_SeChangedFcn ,obj, modelnumber};
    obj.guihandles.(['sigmaOrFactor_' modelnumber]).Position = [0.89 0.33 0.1 0.24];
    obj.guihandles.(['sigmaOrFactor_' modelnumber]).Parent = hParent;
    obj.guihandles.(['useSigmaFactor_' modelnumber]).Parent = obj.guihandles.(['sigmaOrFactor_' modelnumber]);
    obj.guihandles.(['useSigmaFactor_' modelnumber]).Position = [0 32 15 15]; %! the 32 here is used as an flag for the callback 'sigmaOrFactor_SeChangedFcn'
    obj.guihandles.(['useSigma_' modelnumber]).Parent = obj.guihandles.(['sigmaOrFactor_' modelnumber]);
    obj.guihandles.(['useSigma_' modelnumber]).Position = [0 7 15 15];
end
obj.guihandles.(['modelType_' modelnumber]).Enable = 'on';
end

% Callback functions
% Callback for checking a new layer
function layer_callback(a, b, obj, modelNumber)
fitter = obj.fitter;
fitter.resetInit;
fitter.model{modelNumber}.layer = a.Value;
fitter.updateLayer;
fitter.saveInit;
obj.fitter = fitter;
obj.updateLayer;
obj.setPar(['logicalParTableEdited_' obj.name], true); % flag for editing the model settings
end

function modType_callback(a, b, obj, modelNumber)
% Callback for the conversion to an image model.
% Usually only obj.modelType has to be changed unless it is changing
% between image and funciton types.
fitter = obj.fitter;
switch a.String{a.Value}
    case 'image'
        if strcmp(fitter.model{modelNumber}.modelType, 'image')
            warning('The model is already the image type.')
        else
            selection = questdlg('Convert model?','Convert the model to the image type? (Before you convert the model, please make sure you set the right parameters...)','Yes','No','No');
            if strcmp(selection,'Yes')
                prompt = {'Pixel size'};
                dlgtitle = 'Set the pixel size';
                dims = [1 35];
                definput = {'5'};
                pxSize = inputdlg(prompt,dlgtitle,dims,definput);
                pxSize = str2num(pxSize{1});
                % export the image based on the current parameter values and use it
                % as an image model
                
                % save the source model of the image model for later.
                obj.sourceModel{modelNumber} = fitter.model{modelNumber};
                modelImg = fitter.getImage(modelNumber, 'pixelSize', pxSize);
                imgModel = imageModel(modelImg, 'pixelSize', pxSize);
                fitter.changeModel(imgModel,modelNumber);
            end
        end
    otherwise
        if strcmp(fitter.model{modelNumber}.modelType, 'image')
            % do the following if the previous model form is 'image'
            funModel = obj.sourceModel{modelNumber};
            fitter.changeModel(funModel,modelNumber);
            obj.sourceModel{modelNumber} = [];
        else
            fitter.model{modelNumber}.modelType = a.String{a.Value};
        end
end
enableSettings(obj,num2str(modelNumber),fitter);
initmodel(obj, num2str(modelNumber),'skipAddModel',true);
end

% Callback for the model settings edited
function modelSettingsEdited_callback(a,b,obj)
obj.setPar(['logicalParTableEdited_' obj.name], true);
end

% Callback for settingTable
function settingTable_cellEditCallback(a,b,obj, modelnumber)
indices = b.Indices;
data = a.Data;
obj.fitter.linkedGUI = obj;
obj.fitter.loadListener;
for k = 1:obj.fitter.numOfModel
    obj.fitter.model{k}.loadListener;
end
if indices(2)==1
    setting = strtrim(data{indices(1),1});
    value = obj.fitter.getModelInternalSetting(modelnumber,setting);
    data{indices(1),2} = value;
else
    setting = strtrim(data{indices(1),1});
    value = data{indices(1),2};
    obj.fitter.setModelInternalSetting(modelnumber,setting,value);
end
obj.guihandles.(['settingstable_' num2str(modelnumber)]).Data = data;
obj.setPar(['selectedRow_internalSetting_' obj.name '_' num2str(modelnumber)],indices(1));
end

% Callback for settingTable
function settingTable_cellSelectionCallback(a,b,obj, modelnumber)
indices = b.Indices;
if ~isempty(indices)
    obj.setPar(['selectedRow_internalSetting_' obj.name '_' num2str(modelnumber)],indices(1));
end
end

% Callback for internal settings
function addInternalSetting_callback(a,b,obj, modelNumber)
htable = obj.guihandles.(['settingstable_' modelNumber]);
htable.Data = [htable.Data; {[],[]}];
obj.guihandles.(['settingstable_' modelNumber])=htable;
end

function rmInternalSetting_callback(a,b,obj, modelNumber)
htable = obj.guihandles.(['settingstable_' modelNumber]);
data = htable.Data;
selectedRow= obj.getPar(['selectedRow_internalSetting_' obj.name '_' modelNumber]);
try
    data(selectedRow,:) = [];
catch
end
htable.Data = data;
obj.guihandles.(['settingstable_' modelNumber]) = htable;
end

% Action when switching model type
function sigmaOrFactor_SeChangedFcn(a,b,obj, modelNumber)
if b.Source.SelectedObject.Position(2)==32
    obj.fitter.model{str2double(modelNumber)}.fixSigma = true;
    obj.guihandles.(['sigma_fit_' modelNumber]).Enable = 'on';
    obj.guihandles.(['sigmaFactor_fit_' modelNumber]).Enable = 'off';
else
    obj.fitter.model{str2double(modelNumber)}.fixSigma = false;
    obj.guihandles.(['sigma_fit_' modelNumber]).Enable = 'off';
    obj.guihandles.(['sigmaFactor_fit_' modelNumber]).Enable = 'on';
end
end

% Parameters callback
function displayModel_callback(a,b,obj)
% Display the current model
obj.run(obj.getAllParameters, 'onlySetUp',true, 'forceDisplay',true, 'keepParsVal',true);
end

function savePar_callback(a,b,obj,modID)
% Saving the current parameter table
filter = {'*_fitPar.csv'};
[file, path] = uiputfile(filter);
parTable = obj.guihandles.(['partable_' num2str(modID)]).Data;
writecell(parTable, [path file],'Delimiter','tab')
end

function loadPar_callback(a,b,obj,modID)
% Loading and matching a parameter table
filter = {'*_fitPar.csv'};
[file, path] = uigetfile(filter);
importTable = readcell([path file],'Delimiter','\t','LineEnding','\r\n');
dim = size(importTable);
if dim(2) == 8
    importTable = [importTable repmat({char('')},dim(1),1)];
    importTable = cell2table(importTable);
    importTable.Properties.VariableNames = obj.guihandles.(['partable_' num2str(modID)]).ColumnName;
end
refTable = obj.guihandles.(['partable_' num2str(modID)]);
[keyStr,locR,locT] = matchTable(refTable,importTable,{'type','name'},'.');
fig = figure('Name','Import parameter settings');
fig.Position(3:4) = [300 420];
uit = uitable(fig);
uit.Position = [20 40 300-40 420-80];
uit.ColumnName = {'import', 'to'};
ID_tar = strcat(importTable.type,'.',importTable.name);
ID_ref = strcat(refTable.Data(:,6),'.',refTable.Data(:,1));
col2 = cell(size(ID_tar));
col2(locT) = ID_ref(locR);
lEmpty = cellfun(@isempty, col2);
numNotMatched = sum(lEmpty);
if numNotMatched>0
    col2(lEmpty) = cellstr(repmat('none', [numNotMatched 1]));
end
uit.Data = [ID_tar col2];
uit.ColumnFormat = {[],[ID_ref' {'none'}]};
uit.ColumnEditable = [false true];
uit.ColumnWidth = {100, 100};
button = uicontrol(fig,'Style','pushbutton', 'String','Import', 'Callback', {@importPar_callback,obj,modID,uit, importTable, ID_ref});
end

function importPar_callback(a,b,obj,modID,uit, importTable, ID_ref)
refTable = obj.guihandles.(['partable_' num2str(modID)]);
[~,lInput,lOri] = intersect(uit.Data(:,2),ID_ref);
refTable.Data(lOri,[2:5 7:9]) = table2cell(importTable(lInput,[2:5 7:9]));
close(uit.Parent)
end
%

function [keyStr,locR,locT] = matchTable(ref,tar,keys, delimiter)
switch class(ref)
    case 'matlab.ui.control.Table'
        refTab = cell2table(ref.Data);
        refTab.Properties.VariableNames  = ref.ColumnName;
    case 'table'
        refTab = ref;
    otherwise
        error(['ref class "' class(ref) '" not supported. Please use uitable or table instead.'])
        return
end
for k = 1:length(keys)
    if k == 1
        keyStr_ref = refTab.(keys{k});
        keyStr_tar = tar.(keys{k});
    else
        keyStr_ref = strcat(keyStr_ref, delimiter, refTab.(keys{k}));
        keyStr_tar = strcat(keyStr_tar, delimiter, tar.(keys{k}));
        keyStr_ref = deblank(keyStr_ref);
        keyStr_tar = deblank(keyStr_tar);
    end
end
[keyStr,locR,locT] = intersect(keyStr_ref,keyStr_tar);
end

%% # addconverttotab

% ## Define the gui
function pard=guidefconvert(obj)
pard.anchorConvert.object=struct('Style','text','String','');
pard.anchorConvert.position=[10,1];
pard.anchorConvert.Width=2;

pard.addNewRule.object=struct('Style','pushbutton','String','+','Callback',{{@addNewRule_callback,obj,'anchorConvert'}});
pard.addNewRule.position=[10.1,1];
pard.addNewRule.Width=0.2;
pard.addNewRule.Height=0.4;

pard.rmRule.object=struct('Style','pushbutton','String','-','Callback',{{@rmRule_callback,obj,'anchorConvert'}});
pard.rmRule.position=[10.1,1.2];
pard.rmRule.Width=0.2;
pard.rmRule.Height=0.4;

pard.evokeMatchPar.object=struct('Style','pushbutton','String','Match','Callback',{{@evokeParMatcher_callback,obj}});
pard.evokeMatchPar.position=[10.5,4.2];
pard.evokeMatchPar.Width=0.8;
pard.evokeMatchPar.Height=0.8;
end

% ## Match parameters
function hStack = evokeParMatcher_callback(a,b,obj)
obj.setPar('matchPar_fitGUIselected',[]);
fig = figure;
sourceSMLMFit = uicontrol(fig, 'Style', 'popupmenu','Position',[10 380 100 20]);
sourceSMLMFit.String = obj.guihandles.anchorConvert.ColumnFormat{1};

sourceModel = uicontrol(fig, 'Style', 'popupmenu','Position',[10 350 100 20]);
sourceModel.String = {'--'};

sourceSMLMFit.Callback = {@SMLMFitSelect_callback,obj, sourceModel};

matchTable = uitable(fig,'Position',[10 40 250 300]);
allPossibleTargets = obj.guihandles.anchorConvert.ColumnFormat{3};
matchTable.Data = [num2cell(true(size(allPossibleTargets))); allPossibleTargets]';
matchTable.ColumnEditable = [true false];
matchTable.ColumnName = {'Selected', 'Parameter'};

applyMatch = uicontrol(fig, 'Style', 'pushbutton','String','Apply','Position',[10 10 80 30]);
applyMatch.Callback = {@applyMatching_callback,obj, sourceModel,matchTable};

hStack.sourceSMLMFit = sourceSMLMFit;
hStack.sourceModel = sourceModel;
hStack.matchTable = matchTable;
hStack.applyMatch = applyMatch;
end

function SMLMFitSelect_callback(a,b,obj,sourceModel)
selectedFitGUI = a.String{a.Value};
nameAllLoadedEval = obj.locData.SE.processors.eval.guihandles.modules.Data(:,2);
indSelected = ismember(nameAllLoadedEval, selectedFitGUI);
obj.setPar('matchPar_fitGUIselected',indSelected);
numOfModel_selectedFitGUI = obj.locData.SE.processors.eval.processors{indSelected}.fitter.numOfModel;
% export the options
sourceModel.String = cellfun(@num2str, num2cell(1:numOfModel_selectedFitGUI));
end

function parMatcher_callback(a,b,obj)
htable = obj.guihandles.anchorConvert;
data = htable.Data;
selectedRow = obj.getPar(['selectedRowConvert_' obj.name]);
fitter = obj.fitter;
try
    data(selectedRow,:) = [];
catch
end
htable.Data = data;
obj.guihandles.anchorConvert=htable;
end

function applyMatching_callback(a,b,obj, sourceModel,matchTable)
% Match the parameters with the same names.
tabConvert = obj.guihandles.anchorConvert.Data;
indSelected = obj.getPar('matchPar_fitGUIselected');
sourceSMLMModelFit = obj.locData.SE.processors.eval.processors{indSelected}.fitter;
lModel = sourceSMLMModelFit.allParsArg.model == sourceModel.Value;
parType = sourceSMLMModelFit.allParsArg.type(lModel);
parName = sourceSMLMModelFit.allParsArg.name(lModel);
parID_source = join([parType,parName],'.');
lTargetSelected = [matchTable.Data{:,1}];
targetSelected = matchTable.Data(find(lTargetSelected),2);
targetSelected_shourt = regexprep(targetSelected,'m\d\.','');
[lFound,idxSource] = ismember(targetSelected_shourt,parID_source);
ID_target = targetSelected(lFound);
ID_source = parID_source(idxSource(lFound));
prefix = ['pars.m',num2str(sourceModel.Value)];
ID_source = join([cellstr(repmat(prefix,size(ID_source))), ID_source], '.');
fitGUI_source = obj.locData.SE.processors.eval.guihandles.modules.Data{find(indSelected),2};
sizeOfNewRules = size(ID_source);
if isempty(tabConvert)
    tabConvert = [cellstr(repmat(fitGUI_source,sizeOfNewRules)), ID_source, ID_target, cell(sizeOfNewRules)];
else
    tabConvert(end+1:end+sizeOfNewRules(1),1) = cellstr(repmat(fitGUI_source,sizeOfNewRules));
    tabConvert(end-sizeOfNewRules(1)+1:end,2) = ID_source;
    tabConvert(end-sizeOfNewRules(1)+1:end,3) = ID_target;
    tabConvert(end-sizeOfNewRules(1)+1:end,4) = cell(sizeOfNewRules);
end
obj.guihandles.anchorConvert.Data = tabConvert;
end