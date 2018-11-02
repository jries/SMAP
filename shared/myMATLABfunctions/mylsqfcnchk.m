function [allfcns,msg] = lsqfcnchk(funstr,caller,lenVarIn,funValCheck,gradflag)
%LSQFCNCHK Pre- and post-process function expression for FCNCHK.
%   [ALLFCNS,MSG] = LSQFUNCHK(FUNSTR,CALLER,lenVarIn,GRADFLAG) takes
%   the (nonempty) expression FUNSTR from CALLER with LenVarIn extra arguments,
%   parses it according to what CALLER is, then returns a string or inline
%   object in ALLFCNS.  If an error occurs, this message is put in MSG.
%
%   ALLFCNS is a cell array:
%    ALLFCNS{1} contains a flag
%    that says if the objective and gradients are together in one function
%    (calltype=='fungrad') or in two functions (calltype='fun_then_grad')
%    or there is no gradient (calltype=='fun'), etc.
%    ALLFCNS{2} contains the string CALLER.
%    ALLFCNS{3}  contains the objective function
%    ALLFCNS{4}  contains the gradient function (transpose of Jacobian).
%
%    If funValCheck is 'on', then we update the funfcn's (fun/Jacobian) so
%    they are called through CHECKFUN to check for NaN's, Inf's, or complex
%    values. Add a wrapper function, CHECKFUN, to check for NaN/complex
%    values without having to change the calls that look like this:
%    f = funfcn(x,varargin{:}); x is the first argument to CHECKFUN, then
%    the user's function, then the elements of varargin. To accomplish this
%    we need to add the user's function to the beginning of varargin, and
%    change funfcn to be CHECKFUN.
%

%   Copyright 1990-2006 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2009/07/06 20:46:04 $

%    NOTE: we assume FUNSTR is nonempty.
% Initialize
msg='';
allfcns = {};
funfcn = [];
gradfcn = [];

if gradflag
    calltype = 'fungrad';
else
    calltype = 'fun';
end

% {fun}
if isa(funstr, 'cell') && length(funstr)==1
    % take the cellarray apart: we know it is nonempty
    if gradflag
        calltype = 'fungrad';
    end
    [funfcn, idandmsg] = fcnchk(funstr{1},lenVarIn);
    % Insert call to nested function checkfun which calls user funfcn
    if funValCheck
        userfcn = funfcn;
        funfcn = @checkfun; %caller and userfcn are in scope in nested checkfun
    end

    if ~isempty(idandmsg)
        error(idandmsg)
    end

    % {fun,[]}
elseif isa(funstr, 'cell') && length(funstr)==2 && isempty(funstr{2})
    if gradflag
        calltype = 'fungrad';
    end
    [funfcn, idandmsg] = fcnchk(funstr{1},lenVarIn);
    if funValCheck
        userfcn = funfcn;
        funfcn = @checkfun; %caller and userfcn are in scope in nested checkfun
    end
    if ~isempty(idandmsg)
        error(idandmsg)
    end

    % {fun, grad}
elseif isa(funstr, 'cell') && length(funstr)==2 % and ~isempty(funstr{2})

    [funfcn, idandmsg] = fcnchk(funstr{1},lenVarIn);
    if funValCheck
        userfcn = funfcn;
        funfcn = @checkfun; %caller and userfcn are in scope in nested checkfun
    end
    if ~isempty(idandmsg)
        error(idandmsg)
    end
    [gradfcn, idandmsg] = fcnchk(funstr{2},lenVarIn);
    if funValCheck
        userfcn = gradfcn;
        gradfcn = @checkfun; %caller and userfcn are in scope in nested checkfun
    end

    if ~isempty(idandmsg)
        error(idandmsg)
    end
    calltype = 'fun_then_grad';
    if ~gradflag
        warning('optim:lsqfcnchk:IgnoringJacobian', ...
            ['Jacobian function provided but OPTIONS.Jacobian=''off'';\n' ...
            ' ignoring Jacobian function and using finite-differencing.\n' ...
            ' Rerun with OPTIONS.Jacobian=''on'' to use Jacobian function.'])
        calltype = 'fun';
    end

elseif ~isa(funstr, 'cell')  %Not a cell; is a string expression, function name string or inline object
    [funfcn, idandmsg] = fcnchk(funstr,lenVarIn);
    if funValCheck
        userfcn = funfcn;
        funfcn = @checkfun; %caller and userfcn are in scope in nested checkfun
    end

    if ~isempty(idandmsg)
        error(idandmsg)
    end
    if gradflag % gradient and function in one function/M-file
        gradfcn = funfcn; % Do this so graderr will print the correct name
    end
else
    error('optim:lsqfcnchk:InvalidFUN', ...
        ['FUN must be a function or an inline object;\n', ...
        ' or, FUN may be a cell array that contains these type of objects.']);
end

allfcns{1} = calltype;
allfcns{2} = caller;
allfcns{3} = funfcn;
allfcns{4} = gradfcn;
allfcns{5}=[];

%------------------------------------------------------------
    function [f,J] = checkfun(x,varargin)
        % CHECKFUN checks for complex or NaN results from userfcn.
        % Inputs CALLER and USERFCN come from scope of OPTIMFCNCHK.
        % We do not make assumptions about f, or J. For generality, assume
        % they can all be matrices.

        if nargout == 1
            f = userfcn(x,varargin{:});
            if any(any(isnan(f)))
                error('optim:lsqfcnchk:checkfun:NaNFval', ...
                    'User function ''%s'' returned NaN when evaluated;\n %s cannot continue.', ...
                    functiontostring(userfcn),upper(caller));
            elseif ~isreal(f)
                error('optim:lsqfcnchk:checkfun:ComplexFval', ...
                    'User function ''%s'' returned a complex value when evaluated;\n %s cannot continue.', ...
                    functiontostring(userfcn),upper(caller));
            elseif any(any(isinf(f)))
                error('optim:lsqfcnchk:checkfun:InfFval', ...
                    'User function ''%s'' returned Inf or -Inf when evaluated;\n %s cannot continue.', ...
                    functiontostring(userfcn),upper(caller));
            end

        elseif nargout == 2
            [f,J] = userfcn(x,varargin{:});
            if any(any(isnan(f))) || any(any(isnan(J)))
                error('optim:lsqfcnchk:checkfun:NaNFval', ...
                    'User function ''%s'' returned NaN when evaluated;\n %s cannot continue.', ...
                    functiontostring(userfcn),upper(caller));
            elseif ~isreal(f) || ~isreal(J)
                error('optim:lsqfcnchk:checkfun:ComplexFval', ...
                    'User function ''%s'' returned a complex value when evaluated;\n %s cannot continue.', ...
                    functiontostring(userfcn),upper(caller));
            elseif any(any(isinf(f))) || any(any(isinf(J)))
                error('optim:lsqfcnchk:checkfun:InfFval', ...
                    'User function ''%s'' returned Inf or -Inf when evaluated;\n %s cannot continue.', ...
                    functiontostring(userfcn),upper(caller));
            end
        end

    end %checkfun

end % lsqfcnchk