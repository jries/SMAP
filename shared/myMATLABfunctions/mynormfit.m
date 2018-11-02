function [muhat, sigmahat, muci, sigmaci] = normfit(x,alpha,censoring,freq,options)
%NORMFIT Parameter estimates and confidence intervals for normal data.
%   [MUHAT,SIGMAHAT] = NORMFIT(X) returns estimates of the parameters of
%   the normal distribution given the data in X.  MUHAT is an estimate of
%   the mean, and SIGMAHAT is an estimate of the standard deviation.
%
%   [MUHAT,SIGMAHAT,MUCI,SIGMACI] = NORMFIT(X) returns 95% confidence
%   intervals for the parameter estimates.
%
%   [MUHAT,SIGMAHAT,MUCI,SIGMACI] = NORMFIT(X,ALPHA) returns 100(1-ALPHA)
%   percent confidence intervals for the parameter estimates.
%
%   [...] = NORMFIT(X,ALPHA,CENSORING) accepts a boolean vector of the same
%   size as X that is 1 for observations that are right-censored and 0 for
%   observations that are observed exactly.
%
%   [...] = NORMFIT(X,ALPHA,CENSORING,FREQ) accepts a frequency vector of the
%   same size as X.  FREQ typically contains integer frequencies for the
%   corresponding elements in X, but may contain any non-integer
%   non-negative values.
%
%   [...] = NORMFIT(X,ALPHA,CENSORING,FREQ,OPTIONS) specifies control
%   parameters for the iterative algorithm used to compute ML estimates
%   when there is censoring.  This argument can be created by a call to
%   STATSET.  See STATSET('normfit') for parameter names and default values.
%
%   Pass in [] for ALPHA, CENSORING, or FREQ to use their default values.
%
%   With no censoring, SIGMAHAT is computed using the square root of the
%   unbiased estimator of the variance.  With censoring, SIGMAHAT is the
%   maximum likelihood estimate.
%
%   See also NORMCDF, NORMINV, NORMLIKE, NORMPDF, NORMRND, NORMSTAT, MLE, STATSET.

%   References:
%      [1] Evans, M., Hastings, N., and Peacock, B. (1993) Statistical
%          Distributions, 2nd ed., Wiley, 170pp.
%      [2] Lawless, J.F. (1982) Statistical Models and Methods for Lifetime
%          Data, Wiley, New York, 580pp.
%      [3} Meeker, W.Q. and L.A. Escobar (1998) Statistical Methods for
%          Reliability Data, Wiley, New York, 680pp.

%   To compute weighted maximum likelihood estimates (WMLEs) for mu and
%   sigma, you can provide weights, normalized to sum to LENGTH(X), in FREQ
%   instead of frequencies.  In this case, NORMFIT computes the WMLE for
%   mu.  However, when there is no censoring, the estimate computed for
%   sigma is not exactly the WMLE.  To compute the WMLE, multiply the value
%   returned in SIGMAHAT by (SUM(FREQ) - 1)/SUM(FREQ).  This correction is
%   needed because NORMFIT normally computes SIGMAHAT using an unbiased
%   variance estimator when there is no censoring in the data.  When there
%   is censoring, the correction is not needed, since NORMFIT does not use
%   the unbiased variance estimator in that case.

%   Copyright 1993-2008 The MathWorks, Inc.
%   $Revision: 1.1.8.1 $  $Date: 2010/03/16 00:16:28 $

% Illegal data return an error.
if ~isvector(x)
    if nargin < 3
        % Accept matrix data under the 2-arg syntax.  censoring and freq
        % will be scalar zero and one.
        [n,ncols] = size(x); % all columns have same number of data
    else
        error('stats:normfit:InvalidData','X must be a vector.');
    end
else
    n = numel(x); % a scalar -- all columns have same number of data
    ncols = 1;
end

if nargin < 2 || isempty(alpha)
    alpha = 0.05;
end
if nargin < 3 || isempty(censoring)
    censoring = 0; % make this a scalar, will expand when needed
elseif ~isempty(censoring) && ~isequal(size(x), size(censoring))
    error('stats:normfit:InputSizeMismatch',...
          'X and CENSORING must have the same size.');
end
if nargin < 4 || isempty(freq)
    freq = 1; % make this a scalar, will expand when needed
elseif isequal(size(x), size(freq))
    n = sum(freq);
    zerowgts = find(freq == 0);
    if numel(zerowgts) > 0
        x(zerowgts) = [];
        if numel(censoring)==numel(freq), censoring(zerowgts) = []; end
        freq(zerowgts) = [];
    end
else
    error('stats:normfit:InputSizeMismatch',...
          'X and FREQ must have the same size.');
end
if nargin < 5 || isempty(options)
    options = [];
end

classX = class(x);

ncen = sum(freq.*censoring); % a scalar in all cases
nunc = n - ncen; % a scalar in all cases
sumx = sum(freq.*x);

% Weed out cases which cannot really be fit, no data or all censored.  When
% all observations are censored, the likelihood surface is at its maximum
% (zero) for any mu > max(x) at the boundary sigma==0.
if n == 0 || nunc == 0
    muhat = NaN(1,ncols,classX);
    sigmahat = NaN(1,ncols,classX);
    muci = NaN(2,ncols,classX);
    sigmaci = NaN(2,ncols,classX);
    return

% No censoring, find the parameter estimates explicitly.
elseif ncen == 0
    muhat = sumx ./ n;
    if n > 1
        if numel(muhat) == 1 % vector data
            xc = x - muhat;
        else % matrix data
            xc = x - repmat(muhat,[n 1]);
        end
        sigmahat = sqrt(sum(conj(xc).*xc.*freq) ./ (n-1));
    else
        sigmahat = zeros(1,ncols,classX);
    end

    if nargout > 2
        if n > 1
            parmhat = [muhat; sigmahat];
            ci = statnormci(parmhat,[],alpha,x,[],freq);
            muci = ci(:,:,1);
            sigmaci = ci(:,:,2);
        else
            muci = repmat(cast([-Inf; Inf],classX),1,ncols);
            sigmaci = repmat(cast([0; Inf],classX),1,ncols);
        end
    end
    return
end
% Past this point, guaranteed to have only one vector of data, with censoring.

% Not much can be done with Infs, either censored or uncensored.
if ~isfinite(sumx)
    muhat = sumx;
    sigmahat = NaN(classX);
    muci = NaN(2,1,classX);
    sigmaci = NaN(2,1,classX);
    return
end

% When all uncensored observations are equal and greater than all the
% censored observations, the likelihood surface becomes infinite at the
% boundary sigma==0.  Return something reasonable anyway.
xunc = x(censoring==0);
rangexUnc = range(xunc);
if rangexUnc < realmin(classX)
    if xunc(1) == max(x)
        muhat = xunc(1);
        sigmahat = zeros(classX);
        if nunc > 1
            muci = [muhat; muhat];
            sigmaci = zeros(2,1,classX);
        else
            muci = cast([-Inf; Inf], classX);
            sigmaci = cast([0; Inf], classX);
        end
        return
    end
end
% Otherwise the data are ok to fit, go on.

% First, get a rough estimate for parameters using the "least squares" method
% as a starting value...
if rangexUnc > 0
    if numel(freq) == numel(x)
        [p,q] = ecdf(x, 'censoring',censoring, 'frequency',freq);
    else
        [p,q] = ecdf(x, 'censoring',censoring);
    end
    pmid = (p(1:(end-1))+p(2:end)) / 2;
    linefit = polyfit(-sqrt(2)*erfcinv(2*pmid), q(2:end), 1);
    parmhat = linefit([2 1]);

% ...unless there's only one uncensored value.
else
    parmhat = [xunc(1) 1];
end

% The default options include turning statsfminbx's display off.  This
% function gives its own warning/error messages, and the caller can
% turn display on to get the text output from statsfminbx if desired.
options = statset(statset('normfit'), options);
tolBnd = options.TolBnd;
options = optimset(options);
dfltOptions = struct('DerivativeCheck','off', 'HessMult',[], ...
    'HessPattern',ones(2,2), 'PrecondBandWidth',Inf, ...
    'TypicalX',ones(2,1,classX), 'MaxPCGIter',1, 'TolPCG',0.1);

% Maximize the log-likelihood with respect to mu and sigma.
funfcn = {'fungrad' 'normfit' @negloglike [] []};
[parmhat, nll, lagrange, err, output] = ...
         statsfminbx(funfcn, parmhat, [-Inf; tolBnd], [Inf; Inf], ...
                     options, dfltOptions, 1, x, censoring, freq);
if (err == 0)
    % statsfminbx may print its own output text; in any case give something
    % more statistical here, controllable via warning IDs.
    if output.funcCount >= options.MaxFunEvals
        wmsg = 'Maximum likelihood estimation did not converge.  Function evaluation limit exceeded.';
    else
        wmsg = 'Maximum likelihood estimation did not converge.  Iteration limit exceeded.';
    end
    warning('stats:normfit:IterOrEvalLimit',wmsg);
elseif (err < 0)
    error('stats:normfit:NoSolution',...
          'Unable to reach a maximum likelihood solution.');
end

muhat = parmhat(1);
sigmahat = parmhat(2);

if nargout > 2
    parmhat = parmhat(:);
    if numel(freq) == numel(x)
        [nlogL, acov] = normlike(parmhat, x, censoring, freq);
    else
        [nlogL, acov] = normlike(parmhat, x, censoring);
    end
    ci = statnormci(parmhat,acov,alpha,x,censoring,freq);
    muci = ci(:,:,1);
    sigmaci = ci(:,:,2);
end


function [nll,ngrad] = negloglike(parms, x, censoring, freq)
% (Negative) log-likelihood function and gradient for normal distribution.
mu = parms(1);
sigma = parms(2);
cens = (censoring == 1);

z = (x - mu) ./ sigma;
zcen = z(cens);

Scen = .5*erfc(zcen/sqrt(2));
if Scen < realmin(class(x))
    nll = realmax(class(x));
    ngrad = [nll nll];
    return
end
L = -.5.*z.*z - log(sigma);
L(cens) = log(Scen);
nll = -sum(freq .* L);

if nargout > 1
    dlogScen = exp(-.5*zcen.*zcen) ./ (sqrt(2*pi) .* Scen);
    dL1 = z ./ sigma;
    dL1(cens) = dlogScen ./ sigma;
    dL2 = (z.*z - 1) ./ sigma;
    dL2(cens) = zcen .* dlogScen ./ sigma;
    ngrad = [-sum(freq .* dL1) -sum(freq .* dL2)];
end
