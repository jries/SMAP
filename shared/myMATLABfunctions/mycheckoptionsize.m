function checkoptionsize(option, optionsize, numvars, numfunelts)
%CHECKOPTIONSIZE Verify problem size dependent options 
%   CHECKOPTIONSIZE('OPTION', OPTIONSIZE, NUMVARS) verifies that the
%   specified option is of the correct size. Valid values for 'OPTION' are
%   'TypicalX' and 'HessPattern'. OPTIONSIZE is the size of the specified
%   OPTION. NUMVARS specifies the number of free variables (normally
%   numel(X0)).
%
%   CHECKOPTIONSIZE('JacobPattern', OPTIONSIZE, NUMVARS, NUMFUNELTS)
%   verifies that the option 'JacobPattern' is of the correct size.
%   NUMFUNELTS is the number of elements in the vector (or matrix) returned
%   by FUN.
%
%   This is a helper function for the Optimization Toolbox.

%   Copyright 2008 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2008/02/29 13:09:38 $

% Perform size checks
switch lower(option)
    case 'typicalx'
        if prod(optionsize) ~= numvars
            ME = MException('optimlib:checkoptionsize:InvalidSizeOfTypicalX', ...
                'TypicalX must have the same number of elements as X0.');
            throwAsCaller(ME);
        end
    case 'hesspattern'
        expsize = [numvars numvars];
        if ~isequal(optionsize, expsize)
            ME = MException('optimlib:checkoptionsize:InvalidSizeOfHessPattern', ...
                ['User-defined Hessian pattern is not the correct size:\n' ...
                'the matrix HessPattern should be %d-by-%d.'],numvars,numvars);
            throwAsCaller(ME);
        end
    case 'jacobpattern'
        expsize = [numfunelts numvars];
        if ~isequal(optionsize, expsize)
            ME = MException('optimlib:checkoptionsize:InvalidSizeOfJacobPattern', ...
                ['User-defined Jacobian pattern is not the correct size:\n' ...
                'the matrix JacobPattern should be %d-by-%d.'],numfunelts,numvars);
            throwAsCaller(ME);
        end
end

