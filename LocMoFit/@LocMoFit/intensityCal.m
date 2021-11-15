function [totalIntensity, wk, LLctrl] = intensityCal(obj, fitPars, locs, varargin)
% get intensity (likelihood)
%% INTENSITYCAL Get function values (merged models) at given locs
% Call methods 'locsHandler' and 'modelHandler' to get function values at given
% locs
%
% Note:
%   Weight differet models
%   First the weight (weightLayer) will be distributed into every fit layer. Within
%   single layers, the model weights (weightModelLayer) will scale to the locs covered
%   by the models. Each layer has its own offset, summed up to 1, that evenly distributed
%   throughout the ROI. The weight (weightOffset) of offset will be fit.
%
%   Then the final model weight weightModel = (1-weightoffset) * weightLayer *
%   weightModelLayer.
%
%
%

%             This part is for visualizing the
%             obj.allParsArg.value(~obj.allParsArg.fix) = fitPars;
%             fig = figure;
%             viz = axes(fig)
%             obj.plotFreeRot(viz,locs,'lutLocs','redhot','sigma',obj.model{1}.sigma,'pixelSize',obj.model{1}.pixelSize);

%% init
p = inputParser;
p.addParameter('locs2',[]);
p.parse(varargin{:});
locs2 = p.Results.locs2;
fitTheWeights = true; %!!! need to be defined as a changable setting
lPars = {};
mPars = {};
offset = {};
indFit = ~obj.allParsArg.fix;

%% Get intensity values for each model
% Convert and merge the vector fitPars given by the optimizer the fixPars to
% two struct lPars and mPars. Locs are transformed based on lPars and model on
% mPars. Get intensity of the transformed model at the transformed locs.
% Separate fix and fit parameters
% The separation follows the steps bellow:

%%
% # fit lPar
% # remove extracted fitPars
% # fix lPar
% # fit mPar
% # remove extracted fitPars
% # fix mPar

for k = 1:obj.numOfModel
    [lPars,mPars,fitPars] = obj.inputPar2struct(k, fitPars, lPars, mPars, offset);
end
oldLPars = lPars;
lPars = obj.convert2InteralLPar(oldLPars);

for k = 1:obj.numOfModel
    %% handle locs
    % Since all models are aligned to the first model, if the current model
    % is after the first one, then load the already transformed locs.
    if k >1
        locs = obj.getTemp('locsM1');
    end
    
    % fisrt take out the locs that are not in the fit layers
    lLocsInFitLayers = ismember(locs.layer, obj.allModelLayer);
    fn = fieldnames(locs);
    for l = 1:length(fn)
        locs.(fn{l}) = locs.(fn{l})(lLocsInFitLayers,:);
        if ~isempty(locs2)
        end
    end
    
    % transform the locs
    if strcmp(obj.model{k}.modelType,'image')||strcmp(obj.model{k}.modelType,'continuous')
        % this is for models in the image form
        onlyLocpre = false;
    else
        % this is for models in the functional form
        onlyLocpre = true;
    end
    
    newLocs = obj.locsHandler(locs, lPars{k}, k, 'onlyLocpre', onlyLocpre);
    
    if k == 1
        obj.setTemp('locsM1', newLocs);
    end
    % get the size of the largest vector(s) among locs related info
    fn = fieldnames(newLocs);
    numOfFn = length(fn);
    allSize = zeros(numOfFn,2);
    for l = 1:numOfFn
        allSize(l,:) = size(newLocs.(fn{l}));
    end
    maxSize = max(allSize,[],1);
    % Get intensity of one model
    if k == 1
        totalIntensity = zeros(maxSize);
    end
    
    %% Run model handler
    sc = obj.getTemp('currentStep_sigmaCascade');
    
    lSigmaCascade = any(strcmp(obj.model{k}.modelType, {'discrete', 'continuous'}));
    lSigmaCascade = lSigmaCascade||(strcmp(obj.model{k}.modelType, 'image')&&isprop(obj.model{k}, 'linkedLocsModel'));
    
    if lSigmaCascade
        sigF_ori = obj.model{k}.sigmaFactor(1);
%         obj.model{k}.sigmaFactor(1) = sigF_ori*obj.sigmaCascade(1,sc);
        if size(obj.sigmaCascade,1)==2
            obj.model{k}.sigmaFactor(2) = obj.sigmaCascade(2,sc);
        end
    end
    if isempty(mPars)||k>length(mPars)
        if strcmp(obj.model{k}.modelType,'image')
            if isprop(obj.model{k}, 'linkedLocsModel')
                obj.model{k}.checkImg();
            end
            oneIntensity{k} = obj.model{k}.modelHandler(newLocs, []);
        else
            oneIntensity{k} = obj.model{k}.modelHandler(newLocs, [], 'lPars', lPars);
        end
    else
        if strcmp(obj.model{k}.modelType,'image')
            oneIntensity{k} = obj.model{k}.modelHandler(newLocs, mPars{k});
        else
            oneIntensity{k} = obj.model{k}.modelHandler(newLocs, mPars{k}, 'lPars', lPars);
        end
    end
    
    if (~strcmp(obj.model{k}.modelType,'image'))&&lSigmaCascade
        obj.model{k}.sigmaFactor(1) = sigF_ori;
    end
end

%% model weighting
if 1
    
    % This is layer-dependent.
    for layer = 1:obj.numOfLayer
        modLayer = obj.allModelLayer(layer);
        vistedCount = 0;
        for k = obj.numOfModel:-1:1
%             % oneN is a matrix with the same size as the input locs
%             % (number of locs * set of parameters in the optimzation)
%             % remove the locs in the layers that are not used
            if obj.model{k}.layer == modLayer
%                 if vistedCount==0
%                     allN = 0;
%                 end
%                 % Get the count of the locs that are in the model with an intensity value > a certain cutoff.
%                 % count the appearance of the locs in different models
%                 oneN{k} = oneIntensity{k}(ismember(locs.layer, obj.allModelLayer),:)>obj.model{k}.eps; % !!!this should be the value of the normalized gaussian function at sigma of 2
%                 allN = allN + oneN{k};
                vistedCount = vistedCount+1;
                modelInOneLayer{layer}(vistedCount) = k;
            end
        end
        
        if ~fitTheWeights
            % calculate the weighting
            bgN = allN == 0;                    % The number of locs that are purely background
            allN(bgN) = allN(bgN)+1;            % Count the locs that are not in any model as 1 also
            
            for l = length(modelInOneLayer{layer}):-1:1
                k = modelInOneLayer{layer}(l);
                % compensate the model weight
                oneN{k} = sum(oneN{k}./allN,1); % if a loc appear in k models, then the contribution of the locs will be only 1/k
            end
            
            % get the weights
            for l = length(modelInOneLayer{layer}):-1:1
                k = modelInOneLayer{layer}(l);
                wk{k} = oneN{k}./sum(~bgN,1);
            end
        else
            % fit the weighting
            firstVisit = true;
            for l = length(modelInOneLayer{layer}):-1:1
                k = modelInOneLayer{layer}(l);
                if firstVisit
                    allW{layer} = lPars{k}.weight;
                    firstVisit = false;
                    firstInTheLayer=k;
                else
                    %                                 allW{layer} = allW{layer}+lPars{k}.weight-lPars{firstInTheLayer}.weight;
                    allW{layer} = allW{layer}+lPars{k}.weight;
                end
            end
            
            % wk is the weight of models for the same layer
            for l = length(modelInOneLayer{layer}):-1:1
                k = modelInOneLayer{layer}(l);
                if k == firstInTheLayer
                    wk{k} = lPars{k}.weight./allW{layer};
                else
                    %                                 wk{k} = (lPars{k}.weight-lPars{firstInTheLayer}.weight)./allW{layer};
                    wk{k} = (lPars{k}.weight)./allW{layer};
                end
            end            
        end
    end
    
end
if obj.model{1}.dimension == 3
%     offsetUnit = obj.roiSize^-3;
    offsetUnit = 1/(pi*(obj.roiSize/2)^2*obj.roiSize);
else
%     offsetUnit = obj.roiSize^-2;
    offsetUnit = 1/(pi*(obj.roiSize/2)^2);
end
for k = obj.numOfModel:-1:1
    switch 2
        case 1
            totalIntensity = totalIntensity+oneIntensity{k};
            totalIntensity(locs.layer==obj.model{k}.layer,:) = (totalIntensity(locs.layer==obj.model{k}.layer,:)).*oneN{k}./size(allN,1);
        case 2
            % the sum of all the kth weight wk should be 1
            % the offset weight should be between 0-1
            oneIntensity{k}(locs.layer==obj.model{k}.layer,:) = (oneIntensity{k}(locs.layer==obj.model{k}.layer,:).*wk{k}).^obj.model{k}.weight;
            totalIntensity = totalIntensity+oneIntensity{k};
    end
end

% add offset
oldOffset = obj.exportOffset(fitPars);
offset = obj.convert2InteralOffset(oldOffset);

for ch = 1:obj.numOfLayer
    %                 totalIntensity(locs.layer==layer,:) = (((1-offset{90+layer}.weight).*totalIntensity(locs.layer==layer,:) + offset{90+layer}.weight.*offsetUnit))*obj.weightLayer(layer);
%     totalIntensity(locs.layer==layer,:) = (((1-offset{90+layer}.weight).*totalIntensity(locs.layer==layer,:) + (offset{90+layer}.weight.*offsetUnit).^obj.model{k}.weight))*0.5;
    totalIntensity(locs.layer==ch,:) = (((1-offset{90+ch}.weight).*totalIntensity(locs.layer==ch,:) + (offset{90+ch}.weight.*offsetUnit)))*obj.weightLayer(ch);
    % calculate the control
    LLctrl(ch) = log(offsetUnit.*obj.weightLayer(ch))*sum(locs.layer==ch);
end
LLctrl = sum(LLctrl)/sum(ismember(locs.layer,obj.allModelLayer));
% channel not used
used = sum(totalIntensity,2)~=0;
totalIntensity = totalIntensity(used,:);
totalIntensity = double(totalIntensity');

end
