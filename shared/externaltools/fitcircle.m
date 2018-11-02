function [z, r, residual] = fitcircle(x, varargin)
%FITCIRCLE    least squares circle fit
%   
%   [Z, R] = FITCIRCLE(X) fits a circle to the N points in X minimising
%   geometric error (sum of squared distances from the points to the fitted
%   circle) using nonlinear least squares (Gauss Newton)
%       Input
%           X : 2xN array of N 2D points, with N >= 3
%       Output
%           Z : center of the fitted circle
%           R : radius of the fitted circle
%
%   [Z, R] = FITCIRCLE(X, 'linear') fits a circle using linear least
%   squares minimising the algebraic error (residual from fitting system
%   of the form ax'x + b'x + c = 0)
%
%   [Z, R] = FITCIRCLE(X, Property, Value, ...) allows parameters to be
%   passed to the internal Gauss Newton method. Property names can be
%   supplied as any unambiguous contraction of the property name and are 
%   case insensitive, e.g. FITCIRCLE(X, 't', 1e-4) is equivalent to
%   FITCIRCLE(X, 'tol', 1e-4). Valid properties are:
%
%       Property:                 Value:
%       --------------------------------
%       maxits                    positive integer, default 100
%           Sets the maximum number of iterations of the Gauss Newton
%           method
%
%       tol                       positive constant, default 1e-5
%           Gauss Newton converges when the relative change in the solution
%           is less than tol
%
%   [X, R, RES] = fitcircle(...) returns the 2 norm of the residual from 
%   the least squares fit
%
%   Example:
%       x = [1 2 5 7 9 3; 7 6 8 7 5 7];
%       % Get linear least squares fit
%       [zl, rl] = fitcircle(x, 'linear')
%       % Get true best fit
%       [z, r] = fitcircle(x)
%
%   Reference: "Least-squares fitting of circles and ellipses", W. Gander,
%   G. Golub, R. Strebel - BIT Numerical Mathematics, 1994, Springer

% This implementation copyright Richard Brown, 2007, but is freely
% available to copy, use, or modify as long as this line is maintained

error(nargchk(1, 5, nargin, 'struct'))

% Default parameters for Gauss Newton minimisation
params.maxits = 100;
params.tol    = 1e-5;

% Check x and get user supplied parameters
[x, fNonlinear, params] = parseinputs(x, params, varargin{:});

% Convenience variables
m  = size(x, 2);
x1 = x(1, :)';
x2 = x(2, :)';


% 1) Compute best fit w.r.t. algebraic error using linear least squares
% 
% Circle is represented as a matrix quadratic form
%     ax'x + b'x + c = 0
% Linear least squares estimate found by minimising Bu = 0 s.t. norm(u) = 1
%     where u = [a; b; c]

% Form the coefficient matrix
B = [x1.^2 + x2.^2, x1, x2, ones(m, 1)];

% Least squares estimate is right singular vector corresp. to smallest
% singular value of B
[U, S, V] = svd(B);
u = V(:, 4);

% For clarity, set the quadratic form variables
a = u(1);
b = u(2:3);
c = u(4);

% Convert to centre/radius
z = -b / (2*a);
r = sqrt((norm(b)/(2*a))^2 - c/a);

% 2) Nonlinear refinement to miminise geometric error, and compute residual
if fNonlinear
    [z, r, residual] = fitcircle_geometric(x, z, r);
else
    residual = norm(B * u);
end

% END MAIN FUNCTION BODY

% NESTED FUNCTIONS
    function [z, r, residual] = fitcircle_geometric(x, z0, r0)
        % Use a simple Gauss Newton method to minimize the geometric error
        fConverged = false;
        
        % Set initial u
        u     = [z0; r0];
        
        % Delta is the norm of current step, scaled by the norm of u
        delta = inf;
        nIts  = 0;
        
        for nIts = 1:params.maxits
            % Find the function and Jacobian 
            [f, J] = sys(u);
            
            % Solve for the step and update u
            h = -J \ f;
            u = u + h;
            
            % Check for convergence
            delta = norm(h, inf) / norm(u, inf);
            if delta < params.tol
                fConverged = true;
                break
            end
        end
        
        if ~fConverged
            warning('fitcircle:FailureToConverge', ...
                'Gauss Newton iteration failed to converge');
        end
        z = u(1:2);
        r = u(3);
        f = sys(u);
        residual = norm(f);
        
        
        function [f, J] = sys(u)
            %SYS   Nonlinear system to be minimised - the objective
            %function is the distance to each point from the fitted circle
            %contained in u

            % Objective function
            f = (sqrt(sum((repmat(u(1:2), 1, m) - x).^2)) - u(3))';
            
            % Jacobian
            denom = sqrt( (u(1) - x1).^2 + (u(2) - x2).^2 );
            J = [(u(1) - x1) ./ denom, (u(2) - x2) ./ denom, repmat(-1, m, 1)];
        end % sys
        
    end % fitcircle_geometric

% END NESTED FUNCTIONS

end % fitcircle





function [x, fNonlinear, params] = parseinputs(x, params, varargin)
% Make sure x is 2xN where N > 3
if size(x, 2) == 2
    x = x';
end

if size(x, 1) ~= 2
    error('fitcircle:InvalidDimension', ...
        'Input matrix must be two dimensional')
end

if size(x, 2) < 3
    error('fitcircle:InsufficientPoints', ...
        'At least 3 points required to compute fit')
end

% determine whether we are measuring geometric error (nonlinear), or
% algebraic error (linear)
fNonlinear = true;
switch length(varargin)
    % No arguments means a nonlinear least squares with defaul parameters
    case 0
        return
       
    % One argument can only be 'linear', specifying linear least squares
    case 1
        if strncmpi(varargin{1}, 'linear', length(varargin{1}))
            fNonlinear = false;
            return
        else
            error('fitcircle:UnknownOption', 'Unknown Option')
        end
        
    % Otherwise we're left with user supplied parameters for Gauss Newton
    otherwise
        if rem(length(varargin), 2) ~= 0
            error('fitcircle:propertyValueNotPair', ...
                'Additional arguments must take the form of Property/Value pairs');
        end

        % Cell array of valid property names
        properties = {'maxits', 'tol'};

        while length(varargin) ~= 0
            property = varargin{1};
            value    = varargin{2};
            
            % If the property has been supplied in a shortened form, lengthen it
            iProperty = find(strncmpi(property, properties, length(property)));
            if isempty(iProperty)
                error('fitcircle:UnkownProperty', 'Unknown Property');
            elseif length(iProperty) > 1
                error('fitcircle:AmbiguousProperty', ...
                    'Supplied shortened property name is ambiguous');
            end
            
            % Expand property to its full name
            property = properties{iProperty};
            
            switch property
                case 'maxits'
                    if value <= 0
                        error('fitcircle:InvalidMaxits', ...
                            'maxits must be an integer greater than 0')
                    end
                    params.maxits = value;
                case 'tol'
                    if value <= 0
                        error('fitcircle:InvalidTol', ...
                            'tol must be a positive real number')
                    end
                    params.tol = value;
            end % switch property
            varargin(1:2) = [];
        end % while

end % switch length(varargin)

end %parseinputs