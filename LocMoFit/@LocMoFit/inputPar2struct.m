function [lPars, mPars,fitPars] = inputPar2struct(obj, k, fitPars,lPars,mPars,offset)
    % Reshape the input parameters from the optimizer into struct.
    % ---Discription---
    % ---Arguments---
    % Input:
    % obj: SMLMModelFit object.
    % k: the number of the current model
    % fitPars: values of the parameters assigned by the optimizer.
    % Output:
    % lPars: a struct of locs parameters
    % mPars: a struct of model parameters
    % If you would like to perform multi-step fitting, please create one
    % SMLMModelFit object for each step.
    

    indFit = ~obj.allParsArg.fix;
        % get ind
    indModel = obj.allParsArg.model == k;
    indLp = ismember(obj.allParsArg.type,'lPar');
    indMp = ismember(obj.allParsArg.type,'mPar');
    indOs = ismember(obj.allParsArg.type,'offset');

    % fit lPar
    fn = obj.allParsArg.name(indModel&indFit&indLp);
    if ~isempty(fitPars)
        lPars = SMLMModelFit.pars2struct(fn, lPars, fitPars, k,'AddUp',0);
        fitPars = fitPars(:,length(fn)+1:end);
    end
    
    % fix lPar
    fn = obj.allParsArg.name(indModel&~indFit&indLp);
    value = obj.allParsArg.value(indModel&~indFit&indLp);
    lPars = SMLMModelFit.pars2struct(fn, lPars, value', k,'AddUp',0);

    % fit mPar
    fn = obj.allParsArg.name(indModel&indFit&indMp);
    if ~isempty(fitPars)
        mPars = SMLMModelFit.pars2struct(fn, mPars, fitPars, k, 'AddUp', 0);
        fitPars = fitPars(:,length(fn)+1:end);
    end

    % fix mPar
    fn = obj.allParsArg.name(indModel&~indFit&indMp);
    value = obj.allParsArg.value(indModel&~indFit&indMp);
    mPars = SMLMModelFit.pars2struct(fn, mPars, value', k, 'AddUp', 0);
end