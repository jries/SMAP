function stop = optimOutputPar(varargin)
% :func: `optimOutputPar` is used as the outputFcn for optimizers. 
%
% Usage:
%   optimOutputPar(obj,x,optimValues,state)
%   optimOutputPar(obj,optimValues,state)
%
% Args:
%   obj: a LocMoFit object.
%   x: a vector containing the current best parameters.
%   optimValues: a structure containing information of the current optimization iteration.
%   state: the current state of the algorithm, as defined in the optimization toolbox.

obj = varargin{1};
% The different nargin for different optimizer:
% fminsearchbnd passes on 4 elements, while particleswarm 3
if nargin==4
    x = varargin{2};
    optimValues = varargin{3};
    state = varargin{4};
elseif nargin==3
    optimValues = varargin{2};
    state = varargin{3};
    x = optimValues.bestx;
else
    warning('Wrong number of input arguments.')
end
stop = false;
switch state
    case 'init'
        obj.setTemp('optimHistory', obj.parsInit.init(~obj.allParsArg.fix)');
        obj.temp.optimHistory(end+1,:) = x;
    otherwise
        obj.temp.optimHistory(end+1,:) = x;
end
end