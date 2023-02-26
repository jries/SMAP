function [ipoints, residual] = myfindchangepts(x, varargin)
%FINDCHANGEPTS 
%   Reference: 
%   [1] R. Killick, P. Fearnhead, and I.A. Eckley. "Optimal detection of
%       changepoints with a linear computational cost" ArXiv:1101.143v3
%       [stat.ME] 9 Oct 2012
narginchk(1,7);
nargoutchk(0,2);

if nargin > 1
    [varargin{:}] = convertStringsToChars(varargin{:});
end

[statistic, Kmax, Lmin, penalty] = getargs(x, varargin);

wascol = iscolumn(x);
if wascol
  x = x.';
end

if isempty(x) || size(x,2)<Lmin
  cp = zeros(1,0);
  residue = cpnochange(x,statistic);
elseif ~isempty(penalty)
  [cp, residue] = cpmanual(x, statistic, Lmin, penalty);
elseif isempty(Kmax)
  [cp, residue] = cpsingle(x, statistic, Lmin);
  if isempty(cp)
    residue = cpnochange(x,statistic);
  end
else
  [cp, residue] = cpmulti(x, statistic, Lmin, Kmax);
end
  
if nargout==0
  cpplot(x, statistic, cp, residue)
elseif wascol
  ipoints = cp.';
  residual = residue.';
else
  ipoints = cp;
  residual = residue;
end

%-------------------------------------------------------------------------
function [statistic, Kmax, Lmin, penalty] = getargs(x, args)
p = inputParser;
addParameter(p,'Statistic','mean');
addParameter(p,'MaxNumChanges',[]);
addParameter(p,'MinDistance',[]);
addParameter(p,'MinThreshold',[]);

parse(p,args{:});
statistic = p.Results.Statistic;
Kmax = p.Results.MaxNumChanges;
Lmin = p.Results.MinDistance;
penalty = p.Results.MinThreshold;

if ~isempty(penalty) && ~isempty(Kmax)
  error(message('signal:findchangepts:ConflictingParameters', ...
                'MaxNumChanges','MinThreshold'));
elseif ~isempty(Kmax)
  validateattributes(Kmax,{'numeric'},{'integer','nonnegative','scalar'}, ...
                     'findchangepts','MaxNumChanges');
elseif ~isempty(penalty)
  validateattributes(penalty,{'numeric'},{'nonnegative','scalar','nonnan'}, ...
                     'findchangepts','MinThreshold');
end

statistic = validatestring(statistic,{'mean','rms','std','linear'}, ...
                           'findchangepts','Statistic');

if isempty(Lmin)
  if strcmp(statistic,'mean')
    Lmin = 1;
  else
    Lmin = 2;
  end
end

validateattributes(x,{'single','double'},{'real','2d','nonsparse','finite'}, ...
                   'findchangepts','X',1);

validateattributes(Lmin,{'numeric'},{'positive','integer','scalar'}, ...
                   'findchangepts','MinDistance');

if any(strcmp(statistic,{'std','rms','linear'})) && Lmin<2
  error(message('signal:findchangepts:LminTooSmall', ...
                num2str(Lmin),num2str(2),statistic));
end


%-------------------------------------------------------------------------
function residue = cpnochange(x, statistic)
% compute total residual error in the absence of changes
n = size(x,2);
if n==0
  residue = NaN;
elseif strcmp(statistic,'mean')
  residue = n*sum(var(x,1,2));
elseif strcmp(statistic,'rms')
  residue = sum(n*log(sum(x.^2,2)/n));
elseif strcmp(statistic,'std')
  residue = sum(n*log(var(x,1,2)));
elseif strcmp(statistic,'linear')
  residue = sum(n*var(x,1,2) - sum((x-mean(x,2)).*((1:n)-mean(1:n)),2).^2 / (n*var(1:n,1)));
end

%-------------------------------------------------------------------------
function [cp, residue] = cpmanual(x, statistic, Lmin, penalty)
% get upper bound on residual error
residue = cpnochange(x,statistic);

% try one change
[cp, residue1] = cpsingle(x, statistic, Lmin);

% maximum possible penalty
Bmax = residue - residue1;

if isempty(cp) || penalty > Bmax
  % no changes possible
  cp = [];
  return
end

% try the penalty in earnest
[cp, residue] = cppelt(x, penalty, Lmin, statistic);


%-------------------------------------------------------------------------
function [cp, residue] = cpmulti(x, statistic, Lmin, Kmax)

% get upper bound on residual error
resmax = cpnochange(x,statistic);
if Kmax==0
  cp = zeros(1,0);
  residue = resmax;
  return
end

% get lower bound on residual error
[cpmax, minresidue] = cppelt(x, 0, Lmin, statistic);
Pmin = 0;
if Kmax>=numel(cpmax)
  cp = cpmax;
  residue = minresidue;
  return
end


% try one change
[cp, residue] = cpsingle(x, statistic, Lmin);
if isempty(cp)
  residue = cpnochange(x,statistic);
  return
end

% initial penalty
penalty = resmax - residue;
[cp, residue] = cppelt(x, penalty, Lmin, statistic);
Pmax = inf;

% seek lower bound
while numel(cp) < Kmax && residue >= minresidue
  Pmax = penalty;
  resmax = residue;
  cpmax = cp;
  % try reducing the penalty by half
  penalty = 0.5*penalty;
  [cp, residue] = cppelt(x, penalty, Lmin, statistic);
  if numel(cp) > Kmax
    Pmin = penalty;
    minresidue = residue;
  end
end

% seek upper bound
while numel(cp) > Kmax && residue <= resmax && isinf(Pmax)
   Pmin = penalty;
   % try doubling the penalty
   penalty = 2*penalty;
  [cp, residue] = cppelt(x, penalty, Lmin, statistic);
  if numel(cp) < Kmax
    Pmax = penalty;
    cpmax = cp;
    resmax = residue;
  end
end

if numel(cp) == Kmax
  return
end

% search for specified number of changes
penalty = (Pmax + Pmin)/2;
while numel(cp) ~= Kmax && Pmin < penalty && penalty < Pmax
  [cp, residue] = cppelt(x, penalty, Lmin, statistic);
  if numel(cp) < Kmax
    cpmax = cp;
    resmax = residue;
    Pmax = penalty;
  else
    Pmin = penalty;
  end
  penalty = (Pmax + Pmin)/2;
end

% don't exceed Kmax when reporting penalties
if numel(cp) ~= Kmax
  cp = cpmax;
  residue = resmax;
end
