function [xC,Resnorm,FVAL,EXITFLAG,OUTPUT,LAMBDA,JACOB] = ...
    lsqncommon(funfcn,xC,lb,ub,options,defaultopt,caller,initVals, ...
    sizes,flags,mtxmpy,varargin)
%LSQNCOMMON Solves non-linear least squares problems.
%   It contains all the setup code common to both LSQNONLIN and LSQCURVEFIT to 
%   call the different algorithms.

%   Copyright 1990-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2009/08/29 08:31:55 $

lFinite = ~isinf(lb);
uFinite = ~isinf(ub);

[sizes.fRows,sizes.fCols] = size(initVals.F);
sizes.nfun = numel(initVals.F);
initVals.F = initVals.F(:);

diagnostics = isequal(optimget(options,'Diagnostics',defaultopt,'fast'),'on');

largeScaleFlag = strcmp(optimget(options,'LargeScale',defaultopt,'fast'),'on'); % 1 means large-scale, 0 means medium-scale
algorithm = optimget(options,'Algorithm',defaultopt,'fast');
if ~iscell(algorithm)
    initLMparam = 0.01; % Default value
else
    initLMparam = algorithm{2}; % Initial Levenberg-Marquardt parameter
    algorithm = algorithm{1};   % Algorithm string
end

switch algorithm
    case 'trust-region-reflective'
        trustRegion = true;
        lineSearch = false;
    case 'levenberg-marquardt'
        trustRegion = false;
        lineSearch = false;
    case 'lm-line-search'
        % Undocumented Algorithm choice 'lm-line-search'. If it is set, run
        % the old Levenberg-Marquardt line-search code inside nlsq.m
        trustRegion = false;
        lineSearch = true;
        options.LevenbergMarquardt = 'on'; % Needed because it is used in nlsq.m
    otherwise
        error('optim:lsqncommon:InvalidAlgorithm', ...
            ['Invalid choice of option Algorithm for %s. Choose either ''trust-region-reflective''', ...
            ' or ''levenberg-marquardt''.'],upper(caller))
end
levMarq = strcmp(optimget(options,'LevenbergMarquardt',defaultopt,'fast'),'on'); % 0 means Gauss-Newton

% Warn when options conflict. Also warn about options and algorithms being removed.
if trustRegion && ~(largeScaleFlag || levMarq)
    warning('optim:lsqncommon:GNremoval', ...
        ['Gauss-Newton algorithm (selected by LargeScale = ''off'' and ', ...
        'LevenbergMarquardt = ''off'') may be removed in a future release. ', ...
        'Additionally, options LevenbergMarquardt and LargeScale will be ignored ', ...
        'in a future release. Running Gauss-Newton algorithm. Set option Algorithm ', ...
        'to either ''trust-region-reflective'' or ''levenberg-marquardt'' instead.'])
    trustRegion = false;
    lineSearch = true;
elseif trustRegion && ~largeScaleFlag && levMarq
    warning('optim:lsqncommon:AlgorithmConflict', ...
        ['Options LevenbergMarquardt and LargeScale will be ignored in a future release. ',...
        'Running Levenberg-Marquardt algorithm. To run the Levenberg-Marquardt algorithm ',...
        'without this warning, set option Algorithm to ''levenberg-marquardt'' instead.'])
    trustRegion = false;
end

if flags.grad
    % check size of Jacobian
    [Jrows, Jcols] = size(initVals.J);
    if isempty(options.JacobMult) 
        % Not using 'JacobMult' so Jacobian must be correct size
        if Jrows ~= sizes.nfun || Jcols ~= sizes.nVar
            error('optim:lsqncommon:InvalidJacSize',['User-defined Jacobian is not the correct size:\n' ...
                     'the Jacobian matrix should be %d-by-%d.'],sizes.nfun,sizes.nVar)
        end
    end
else
    Jrows = sizes.nfun; 
    Jcols = sizes.nVar;   
end

% trust-region-reflective and not enough equations -- switch to other methods
if trustRegion && sizes.nfun < sizes.nVar
    if (isempty(lb(lFinite)) && isempty(ub(uFinite)))
        warning('optim:lsqncommon:SwitchToLineSearch', ...
            ['Trust-region-reflective algorithm requires at least as many equations' ...
            ' as variables; using Levenberg-Marquardt algorithm instead.'])
    end
    trustRegion = false; % Call levenbergMarquardt
    lineSearch = false;
% Leveberg-Marquardt or line-search with bounds and enough equations,
% switch to trust-region
elseif ~trustRegion && (~isempty(lb(lFinite)) || ~isempty(ub(uFinite))) && sizes.nfun >= sizes.nVar
    warning('optim:lsqncommon:SwitchToLargeScale', ...
        ['Levenberg-Marquardt and Gauss-Newton algorithms do not handle bound' ...
        ' constraints; using trust-region-reflective algorithm instead.'])
    trustRegion = true;
    lineSearch = false;
end
% Can't handle this one:
if (~isempty(lb(lFinite)) || ~isempty(ub(uFinite))) && sizes.nfun < sizes.nVar
    error('optim:lsqncommon:ProblemNotHandled', ...
        ['Levenberg-Marquardt and Gauss-Newton algorithms do not handle bound' ...
        ' constraints and trust-region-reflective algorithm requires at least' ...
        ' as many equations as variables; aborting.']);
end

if diagnostics > 0
    % Do diagnostics on information so far
    constflag=0;gradconstflag=0;non_eq=0;non_ineq=0;lin_eq=0;lin_ineq=0;
    confcn{1}=[]; c=[]; ceq=[]; cGRAD=[]; ceqGRAD=[];
    hessflag = 0; HESS=[];
    if trustRegion
        OUTPUT.algorithm = 'trust-region-reflective';
    elseif lineSearch
        OUTPUT.algorithm = 'medium-scale: line-search';
    else
        OUTPUT.algorithm = 'Levenberg-Marquardt';
    end
    diagnose(caller,OUTPUT,flags.grad,hessflag,constflag,gradconstflag,...
        ~largeScaleFlag,options,defaultopt,xC(:),non_eq,non_ineq,...
        lin_eq,lin_ineq,lb,ub,funfcn,confcn,initVals.F,initVals.J,HESS,c,ceq,cGRAD,ceqGRAD);
end

% Prepare strings to give feedback to users on options they have or have not set.
% These are used in the exit messages.
optionFeedback = mycreateOptionFeedback(options);

% Execute algorithm
if trustRegion
    if ~flags.grad % provide sparsity of Jacobian if not provided.
        Jstr = optimget(options,'JacobPattern',[]);
        if isempty(Jstr)  
            % Put this code separate as it might generate OUT OF MEMORY error
            Jstr = sparse(ones(Jrows,Jcols));
        elseif ischar(Jstr) 
            % options.JacobPattern is the default: 'sparse(ones(jrows,jcols))'
            Jstr = sparse(ones(Jrows,Jcols));
        else 
            % Pattern matrix  - other datatypes (cell-array, struct) are checked in optimset and its
            % helper functions
            checkoptionsize('JacobPattern', size(Jstr), sizes.nVar, sizes.nfun);
        end
    else
        Jstr = [];
    end
    % Set MaxFunEvals appropriately for trust-region-reflective
    defaultopt.MaxFunEvals = '100*numberOfVariables';
    [xC,FVAL,LAMBDA,JACOB,EXITFLAG,OUTPUT,msgData]=...
        mysnls(funfcn,xC,lb,ub,flags.verbosity,options,defaultopt,initVals.F,initVals.J,caller, ...
        Jstr,flags.computeLambda,mtxmpy,flags.detailedExitMsg,optionFeedback,varargin{:});
elseif lineSearch
    % Set MaxFunEvals appropriately for line-search methods
    defaultopt.MaxFunEvals = '100*numberOfVariables';
    [xC,FVAL,JACOB,EXITFLAG,OUTPUT,msgData] = ...
        nlsq(funfcn,xC,flags.verbosity,options,defaultopt,initVals.F,initVals.J,caller,varargin{:});
    LAMBDA.upper = zeros(sizes.nVar,1); 
    LAMBDA.lower = zeros(sizes.nVar,1);
else
    % Set MaxFunEvals appropriately for Levenberg-Marquardt
    defaultopt.MaxFunEvals = '200*numberOfVariables';
    [xC,FVAL,JACOB,EXITFLAG,OUTPUT,msgData] = ...
        levenbergMarquardt(funfcn,xC,flags.verbosity,options,defaultopt,initVals.F,initVals.J, ...
        caller,initLMparam,flags.detailedExitMsg,optionFeedback,varargin{:});
    LAMBDA.upper = zeros(sizes.nVar,1); 
    LAMBDA.lower = zeros(sizes.nVar,1);
end

Resnorm = FVAL'*FVAL;
if ~lineSearch % If nlsq was not called
    OUTPUT.message = createExitMsg(msgData{:});
else % nlsq was called
    OUTPUT.message = msgData;
    if flags.verbosity > 0
        disp(OUTPUT.message);
    end
end

% Reset FVAL to original shapes
FVAL = reshape(FVAL,sizes.fRows,sizes.fCols);

