function fit(obj, locs, varargin)
%% FIT Perform fitting based on the option values
% :func:`fit` performs fitting based on the options.
%
% Usage:
%   fit(obj, locs, varargin)
%
% Args:
%   obj: a SMLMModelFit object.
%   locs: a struct with fields of xnm, ynm, znm, locprecnm.
%   Name-value pairs:
%       *'locs2': a struct with fields of xnm, ynm, znm, locprecnm. This is the
% second set of locs for particle fusion.


%% Deal with the varargin
obj.locs = locs;
p = inputParser;
p.addParameter('locs2', [])
p.addParameter('controlLogLikelihood', obj.getAdvanceSetting('controlLogLikelihood'))
p.addParameter('confidenceInterval', obj.getAdvanceSetting('confidenceInterval'))
p.addParameter('skipFit', false)
p.parse(varargin{:})
p = p.Results;
locs2 = p.locs2;                % the second set of locs
obj.useCompensation = true;

%% Set up the optimizer
% for compiled version
@fminsearchbnd;
@fmincon;
@particleswarm;

solverFun = str2func(obj.solver.SolverName);
if ~isempty(obj.solver.SolverOptions)&&~p.skipFit
    obj.status = 'initial';
    useOptimset = {'fmincon', 'fminsearchbnd'};
    % Change the class of the solver optoins according to the items
    % themsleves.
    ind = ismember(obj.solver.SolverOptions(1:2:end), {'TolFun', 'TolX','MaxIter','MaxFunEvals'});
    idx = find(ind);
    if sum(ind)>0
        for k = length(idx):-1:1
            indSub(k) = ~isnumeric(obj.solver.SolverOptions{idx(k)*2});
        end
        idx = idx(indSub);
        obj.solver.SolverOptions(idx*2) = num2cell(str2double(obj.solver.SolverOptions(idx*2)));
    end
    
    % Deal with the data type of the solverOptions
    %   Convert all that can be converted to numeric
    solverOptions = obj.solver.SolverOptions;
    strInd = find(cellfun(@(x)isstr(x),solverOptions));
    for k = 1:length(strInd)
        val = str2num(solverOptions{strInd(k)});
        if ~isempty(val)
            solverOptions{strInd(k)} = val;
        end
    end
    indF = strcmp(solverOptions, 'false');
    if any(indF)
        solverOptions{indF} = false;
    end
    indT = strcmp(solverOptions, 'true');
    if any(indT)
        solverOptions{indT} = true;
    end
    
    % Deal with the plotFcn
    plotFun = [];
    indPlotFcn = find(strcmp(solverOptions,'PlotFcn'));
    if ~isempty(indPlotFcn)
        if strcmp(solverOptions{indPlotFcn+1},'on')
            plotFun = solverOptions(indPlotFcn:indPlotFcn+1);
        end
        solverOptions(indPlotFcn:indPlotFcn+1) = [];
    end
    outputFcn = [];
    indOutputFcn = find(strcmp(solverOptions,'OutputFcn'));
    if ~isempty(indOutputFcn)
        if strcmp(solverOptions{indOutputFcn+1},'on')
            outputFcn = solverOptions(indOutputFcn:indOutputFcn+1);
        end
        solverOptions(indOutputFcn:indOutputFcn+1) = [];
    end
    
    if ismember(obj.solver.SolverName, useOptimset)
        solverOption = optimset(solverOptions{:});
    else
        solverOption = optimoptions(obj.solver.SolverName, solverOptions{:});
    end
    %     solverOption.PlotFcns = {@optimplotx};
    
    % Deal with the plotFcn/outputFcn
    if ~isempty(plotFun)
        solverOption.PlotFcns = @obj.optimoPlotMod;
    end
    if ~isempty(outputFcn)
        solverOption.OutputFcn = @obj.optimOutputPar;
    end
else
    solverOption = {};
end

%% Set the init values of the parameters based on the rules
% The rules are defined by the users to convert other information to the initial
% guess, and upper/lower bounds etc.
if ~p.skipFit
    if ~isempty(obj.converterRules)                                                     % if the converter rules exist
        obj.convertNow(locs);
    end
end

if ~p.skipFit
    obj.preFittingConversion;
end

for sc = 1:size(obj.sigmaCascade,2) % This is for the sigma cascading
    obj.currentCascadeStep = sc;
    [lb, ub, init, minVal, maxVal] = obj.prepFit;                                   % get settings for fit parameters
    if sc > 1
        init = parBestFit';
    end
    obj.setTemp('currentStep_sigmaCascade', sc);                                    % letting the intensityCal know the current step
    lb = lb*obj.sigmaCascade(1,sc);
    ub = ub*obj.sigmaCascade(1,sc);
    % use the current values as init
    obj.parsInit.init = obj.allParsArg.value;
    % parameters of the best fit
    finalLb = lb'+init';                                                            % the lb and ub are relative lower and upper bounds that centered at the initial guess init
    finalUb = ub'+init';
    
    % Check the bounds are out of meaningful range or not.
    % The meaningful range is defined by user in the model functions.
    % If the bounds are out of the range, then reset the bounds to the min
    %   and max meaningful values
    lSet2Min = finalLb<minVal';
    lSet2Max = finalUb>maxVal';
    finalLb(lSet2Min) = minVal(lSet2Min);
    finalUb(lSet2Max) = maxVal(lSet2Max);
    
    %% Get locs info
    obj.getLocsInfo
    
    %% Get compensationFactor
    % This factor compensate the number of locs between different channels
    switch obj.getAdvanceSetting('layerNorm')
        case 'on'
            compensationFactor = zeros(size(locs.layer))';
            for k = 1:obj.numOfLayer
                compensationFactor(locs.layer==k) = obj.compensationFactor(k);
            end
        case 'off'
            for k = 1:obj.numOfLayer
                compensationFactor(locs.layer==k) = 1;
            end
    end
    compensationFactor(~ismember(locs.layer, obj.allModelLayer)) = [];
    %% Run the optimization
    
    switch obj.objFunType
        case 'likelihood'
            objFun = @(fitPars)-obj.loglikelihoodFun(fitPars, compensationFactor, locs);
        case 'correlation'
            objFun = @(fitPars)-sum(obj.intensityCal(fitPars,locs),2);
    end
    
    indFit = ~obj.allParsArg.fix;
    if ~p.skipFit
        
        switch obj.getAdvanceSetting('runtime')
            case 'on'
                tic
            case 'off'
%               'Do nothing.'
        end
        switch obj.solver.SolverName
            case 'fminsearchbnd'
                [parBestFit,mLLfit] = solverFun(@(fitPars)objFun(fitPars),...
                    init',...
                    finalLb,...
                    finalUb,...
                    solverOption);
            case 'particleswarm'
                solverOption.InitialSwarmMatrix = init';
                [parBestFit,mLLfit] = solverFun(@(fitPars)objFun(fitPars),...
                    length(lb),...
                    finalLb,...
                    finalUb,...
                    solverOption);
            case 'fmincon'
                [parBestFit,mLLfit] = solverFun(@(fitPars)objFun(fitPars),...
                    init',...
                    [],...
                    [],...
                    [],...
                    [],...
                    finalLb,...
                    finalUb,...
                    [],...
                    solverOption);
        end
        switch obj.getAdvanceSetting('runtime')
            case 'on'
                runtime = toc;
            case 'off'
%               'Do nothing.'
        end
        obj.postFittingConversion(parBestFit);
    else
        oldFitInfo = obj.fitInfo;
        parBestFit = obj.allParsArg.value(indFit)';
        mLLfit = -obj.fitInfo.LLfit*sum(ismember(locs.layer,obj.allModelLayer));
        if isfield(obj.fitInfo, 'runtime')
            runtime = obj.fitInfo.runtime;
        end
    end
end
obj.currentCascadeStep = 1; % reset the currenct step
obj.allParsArg.value(indFit) = parBestFit;            % only update the fit parameters

if ~p.skipFit
    %% Compensation of the pixel-size-depedent offsets
    switch obj.model{1}.modelType
        % Since all models are relative to model 1, only its xyz need to be
        % corrected.
        case {'continuous','image'}
            if obj.dataDim==2
                parName = {'x','y'};
            else
                parName = {'x','y','z'};
            end
            
            for k = 1:length(parName)
                [val,ind] = obj.getVariable(['m1.' parName{k}]);
                if ~obj.allParsArg.fix(ind)
                    obj.setParArg(['m1.lPar.' parName{k}], 'value', val+obj.model{1}.pixelSize);
                end
            end
        otherwise
            % do nothing
    end
    obj.status = 'finished';
end

%% Calculate log-likelihood values
% save the log liklihood and weighting factors
% fitInfo is created here
if ~p.skipFit
    obj.fitInfo = [];
    numFittedLocs = sum(ismember(locs.layer,obj.allModelLayer));
    obj.fitInfo.LLfit = -mLLfit/numFittedLocs;
    obj.fitInfo.numOfLocsPerLayer = obj.numOfLocsPerLayer;
    obj.fitInfo.BGDensity = obj.getBGDensity;
    obj.fitInfo.AIC = 2*obj.numOfFreePar+2*mLLfit;
    obj.fitInfo.AICc = obj.fitInfo.AIC+2*(obj.numOfFreePar^2+obj.numOfFreePar)/(numFittedLocs-obj.numOfFreePar-1);
    obj.fitInfo.normAICc = obj.fitInfo.AICc/numFittedLocs;
    if isfield(obj.temp, 'optimHistory')
        obj.fitInfo.optimHistory = obj.getTemp('optimHistory');
        obj.rmTemp('optimHistory');
    end
    if exist('runtime','var')
        obj.fitInfo.runtime = runtime;
    end
    %% Calculate control log-likelihood values
    switch p.controlLogLikelihood
        case 'none'
            obj.fitInfo.LLExpDist = [];
        case 'expected'
    %         obj.fitInfo.LLExp = obj.getELL(parBestFit,compensationFactor,5);
            obj.fitInfo.LLExpDist = obj.getLLExpDist(1000);
            obj.fitInfo.LLExp = mean(obj.fitInfo.LLExpDist);
            obj.fitInfo.LLExpStd = std(obj.fitInfo.LLExpDist);
            obj.fitInfo.LLZScore = (obj.fitInfo.LLfit-obj.fitInfo.LLExp)./obj.fitInfo.LLExpStd;
        case 'overfitted'
            obj.fitInfo.LLOF = obj.getOFLL(compensationFactor);
        case 'both'
            obj.fitInfo.LLExp = obj.getELL(parBestFit,compensationFactor,2);
            obj.fitInfo.LLOF = obj.getOFLL(compensationFactor);
    end
else
    obj.fitInfo = oldFitInfo;
end

%% Variations of parameters
try
    if strcmp(p.confidenceInterval, 'on')
        %% Numerical Hessian matrix
        % This is for a rough estimation of the parameter ranges
        H = hessian(@(fitPars)objFun(fitPars), parBestFit);
        invH = inv(H);
        par_std_int = abs(sqrt(diag(invH)));
        
        %% Approximate the log-likelihood function with a quadratic form
        % Here the quadratic form is defined as LL = x'Hx+a.
        % H is the hessian matrix that we want to estimate. x is a
        % np-element vector of parameter values (centered to the best
        % parameters), where 'np' is the number of fitted parameters. LL is
        % a scalar of the corresponding log-likelihood value.
        
        numOfPar = length(parBestFit);
        iter = 0;
        par_std = 1i; % set par_std as a imaginary number
        
        % Fit the quadratic form to [numOfSampling] points around the best
        % parameters.
        while iter<10&&~isreal(par_std)
            numOfSampling = 1000;
            par = parBestFit+...
                rand(numOfSampling,length(parBestFit)).*...
                2.*par_std_int'-par_std_int';
            allLL = [];
            
            for k_LL = size(par,1):-1:1
                allLL(k_LL) = objFun(par(k_LL,:));
                while ~isreal(allLL(k_LL))
                    % repeat until no imaginary LL caseud by out-of-domain
                    % sampling.
                    par(k_LL,:) = parBestFit+rand(1,length(parBestFit)).*...
                        2.*par_std_int'-par_std_int';
                    allLL(k_LL) = objFun(par(k_LL,:));
                end
            end
            
            % Fitting
            quaTerms = getQuadraticTerms(numOfPar);
            model = polyfitn(par-parBestFit,allLL,quaTerms);
            H = zeros(numOfPar);
            for k_term = 1:length(quaTerms)
                ind = find(quaTerms(k_term,:));
                if length(ind)==1
                    H(ind,ind) = model.Coefficients(k_term);
                elseif length(ind)==2
                    H(ind(1),ind(2)) = model.Coefficients(k_term)./2;
                    H(ind(2),ind(1)) = model.Coefficients(k_term)./2;
                end
            end
            
            
            if false
                % Model vs data plot for a visual inspection
                for k_sampling = 1:numOfSampling
                    allLL_qf(k_sampling) = (par(k_sampling,:)-parBestFit)*H*(par(k_sampling,:)-parBestFit)';
                end
                k = 8; figure; plot(par(:,k),allLL_qf+model.Coefficients(end), ' .'); hold on; plot(par(:,k),allLL, ' .'); hold off
            end
            
            % Calculations of the parameter stds. Quite often the sampling
            % windowis are too small to have a convex parameter-to-LL
            % function. In this case, the windows are expanded [winExpand]
            % times.
            winExpand = 1.5;
            invH = inv(H);
            indBad = diag(invH)<0;
            par_std_int(indBad) = par_std_int(indBad).*winExpand;
            par_std = sqrt(diag(invH));
            iter = iter+1;
        end
        
        % Output results
        parID = obj.getAllParId();
        CI_ub = parBestFit+1.96*par_std';
        CI_lb = parBestFit-1.96*par_std';
        obj.fitInfo.CI.parameter = parID(indFit)';
        obj.fitInfo.CI.interval = [CI_ub; CI_lb];
        obj.fitInfo.CI.std = par_std';
        obj.fitInfo.CI.adjR2 = model.AdjustedR2;
    end
catch ME
    warning(['Site ' num2str(obj.linkedGUI.site.indList) ': cannot estimate the hessian matrix properly.'])
    disp(getReport(ME, 'extended', 'hyperlinks', 'on'))
    obj.fitInfo.CI = [];
end
for k=obj.numOfLayer:-1:1
    obj.fitInfo.numOfLocsPerLayer_BGFree(k) = obj.fitInfo.numOfLocsPerLayer(k) * (1-obj.getVariable(['pars.m9' num2str(k) '.offset.weight']));
end
[~,obj.fitInfo.weightModel,obj.fitInfo.LLctrl] = obj.intensityCal(parBestFit, locs);
for k = 1:obj.numOfModel
    if ~strcmp(obj.model{k}.modelType, 'image')
        [~, modelPar_internal] = obj.model{k}.getPoint(obj.exportPars(k,'mPar'));
        obj.fitInfo.modelPar_internal{k} = modelPar_internal;
    end
end
end

function quaTerms = getQuadraticTerms(numOfPar)
    allTerms = permn([1 2 0], numOfPar);
    ind = sum(allTerms,2)==2;
    allTerms = allTerms(ind,:);
    quaTerms = [allTerms; zeros(1,size(allTerms,2))];
end