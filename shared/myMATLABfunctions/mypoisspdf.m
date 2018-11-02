function y = poisspdf(x,lambda)
%POISSPDF Poisson probability density function.
%   Y = POISSPDF(X,LAMBDA) returns the Poisson probability density 
%   function with parameter LAMBDA at the values in X.
%
%   The size of Y is the common size of X and LAMBDA. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%
%   Note that the density function is zero unless X is an integer.
%
%   See also POISSCDF, POISSFIT, POISSINV, POISSRND, POISSTAT, PDF.

%   References:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.1.22. 
%      [2]  C. Loader, "Fast and Accurate Calculations of Binomial
%      Probabilities", 2000.

%   Copyright 1993-2010 The MathWorks, Inc.
%   $Revision: 1.1.8.3 $  $Date: 2010/05/10 17:59:15 $

if nargin <  2, 
    error('stats:poisspdf:TooFewInputs','Requires two input arguments.'); 
end

[errorcode x lambda] = mydistchck(2,x,lambda);

if errorcode > 0
    error('stats:poisspdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

if isa(x,'single') || isa(lambda,'single')
   y = zeros(size(x),'single');
else
   y = zeros(size(x));
end

if ~isfloat(x)
   x = double(x);
end

y(lambda < 0) = NaN;
y(isnan(x) | isnan(lambda)) = NaN;
y(x==0 & lambda==0) = 1;

k = find(x >= 0 & x == round(x) & lambda > 0);

if ~isempty(k)
    x = x(k);
    lambda = lambda(k);
 
    smallx = x<=lambda*realmin;
    y(k(smallx)) = exp(-lambda(smallx));
    
    largex = lambda<x*realmin;
    y(k(largex)) = exp(-lambda(largex) + x(largex).*log(lambda(largex)) ...
        - gammaln(x(largex)+1));
    
    other = ~smallx & ~largex;
    lnsr2pi = 0.9189385332046727; % log(sqrt(2*pi))
    y(k(other)) = exp(-lnsr2pi -0.5*log(x(other)) - stirlerr(x(other)) ...
        - binodeviance(x(other),lambda(other)));
end
