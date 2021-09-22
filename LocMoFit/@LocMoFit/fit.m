function fit(obj, locs, varargin)
%% FIT Perform fitting based on the option values
% :func:`fit` performs fitting based on the options
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
p.addParameter('skipFit', false)
p.parse(varargin{:})
p = p.Results;
locs2 = p.locs2;                % the second set of locs
obj.useCompensation = true;

%% Set up the optimizer
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
    
    solverOptions = obj.solver.SolverOptions;
    
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

for sc = 1:length(obj.sigmaCascade) % This is for the sigma cascading
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
    finalLb(lSet2Min)=minVal(lSet2Min);
    finalUb(lSet2Max)=maxVal(lSet2Max);
    
    %% Get locs info
    obj.getLocsInfo
    
    %% Get compensationFactor
    % This factor compensate the number of locs between different channels
    compensationFactor = zeros(size(locs.layer))';
    for k = 1:obj.numOfLayer
        compensationFactor(locs.layer==k) = obj.compensationFactor(k);
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
        switch obj.solver.SolverName
            case 'fminsearchbnd'
                [parBestFit,mLLfit] = solverFun(@(fitPars)objFun(fitPars),...
                    init',...
                    finalLb,...
                    finalUb,...
                    solverOption);
            case 'particleswarm'
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
    else
        parBestFit = obj.allParsArg.value(indFit)';
        mLLfit = -obj.fitInfo.LLfit*sum(ismember(locs.layer,obj.allModelLayer));
    end
end
obj.currentCascadeStep = 1; % resent the currenct step
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

%% Calculate control log-likelihood values

switch p.controlLogLikelihood
    case 'none'
    case 'expected'
        obj.fitInfo.LLExp = obj.getELL(parBestFit,compensationFactor,5);
        obj.fitInfo.LLExpDist = obj.getLLExpDist(1000);
        obj.fitInfo.LLExpMean = mean(obj.fitInfo.LLExpDist);
        obj.fitInfo.LLExpStd = std(obj.fitInfo.LLExpDist);
    case 'overfitted'
        obj.fitInfo.LLOF = obj.getOFLL(compensationFactor);
    case 'both'
        obj.fitInfo.LLExp = obj.getELL(parBestFit,compensationFactor,2);
        obj.fitInfo.LLOF = obj.getOFLL(compensationFactor);
    otherwise
        obj.fitInfo.LLExpDist = [];
end


if 0
    % hessian related values
    parID = obj.getAllParId();
    par_std = sqrt(hessdiag(@(fitPars)objFun(fitPars), parBestFit));
    CI_ub = parBestFit+1.96*par_std;
    CI_lb = parBestFit-1.96*par_std;
    obj.fitInfo.CI.parameter = parID(indFit)';
    obj.fitInfo.CI.interval = [CI_ub; CI_lb];
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