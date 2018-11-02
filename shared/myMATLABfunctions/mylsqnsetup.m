function [funfcn,mtxmpy,flags,sizes,funValCheck,xstart,lb,ub,EXITFLAG,Resnorm,FVAL, ...
    LAMBDA,JACOB,OUTPUT,earlyTermination] = lsqnsetup(FUN,xC,LB,UB,options, ...
    defaultopt,caller,callerNargout,lengthVarargin)
%LSQNSETUP Helper function for lsqcurvefit and lsqnonlin.
% This function performs setup tasks common to both solvers.

%   Copyright 2008-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2009/07/06 20:46:07 $

display = optimget(options,'Display',defaultopt,'fast');
flags.detailedExitMsg = ~isempty(strfind(display,'detailed'));
switch display
case {'off','none'}
    flags.verbosity = 0;
case {'iter','iter-detailed'}
    flags.verbosity = 2;
case {'final','final-detailed'}
    flags.verbosity = 1;
case 'testing'
    flags.verbosity = Inf;
otherwise
    flags.verbosity = 1;
end

funValCheck = strcmp(optimget(options,'FunValCheck',defaultopt,'fast'),'on');
flags.grad =  strcmp(optimget(options,'Jacobian',defaultopt,'fast'),'on');

if callerNargout > 5
   flags.computeLambda = true;
else 
   flags.computeLambda = false;
end

% Process objective function
if ~isempty(FUN)  % will detect empty string, empty matrix, empty cell array
    funfcn = mylsqfcnchk(FUN,caller,lengthVarargin,funValCheck,flags.grad);
else
    mexcptn = MException(['optim:' caller ':InvalidFUN'], ...
        'FUN must be a function handle; or, FUN may be a cell array that contains function handles.');
    throwAsCaller(mexcptn);
end

% Process function JacobMult
mtxmpy = optimget(options,'JacobMult',defaultopt,'fast'); 
% Check if name clash
myfunctionNameClashCheck('JacobMult',mtxmpy,'atamult','optim:lsqnsetup:JacobMultNameClash');

% Use internal Jacobian-multiply function if user does not provide JacobMult function 
% or options.Jacobian is off
if isempty(mtxmpy) || (~strcmpi(funfcn{1},'fungrad') && ~strcmpi(funfcn{1},'fun_then_grad'))
    mtxmpy = @atamult;
end

xstart = xC(:); % operate on column-vector xstart
sizes.nVar = numel(xstart);
[xstart,lb,ub,msg] = mycheckbounds(xstart,LB,UB,sizes.nVar);
if ~isempty(msg)
    EXITFLAG = -2;
    [Resnorm,FVAL,LAMBDA,JACOB] = deal([]);
    OUTPUT.firstorderopt = [];
    OUTPUT.iterations = 0;
    OUTPUT.funcCount = 0;
    OUTPUT.cgiterations = [];
    OUTPUT.algorithm = ''; % Not known at this stage
    OUTPUT.message = msg;
    if flags.verbosity > 0
        disp(msg)
    end
    earlyTermination = true; % exit this helper function and return
else
    % All outputs need to be assigned to a value
    Resnorm = []; FVAL = []; LAMBDA = []; JACOB = []; OUTPUT = [];
    EXITFLAG = []; earlyTermination = false;
    
    % If components of initial x not within bounds, set those components
    % of initial point to a "box-centered" point
    xinitOutOfBounds_idx = xstart < lb | xstart > ub;
    if any(xinitOutOfBounds_idx)
        xstart = startx(ub,lb,xstart,xinitOutOfBounds_idx);
    end
end


