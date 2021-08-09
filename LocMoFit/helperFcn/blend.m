function foo = blend(f1,f2,location,distance)
% BLEND  Uses hyperbolic tangent (tanh) to smoothly blend two functions or
% smoothly combine a piecewise function into a single function handle.
% 
%   F = BLEND(f1,f2,LOC,DIST) returns function handle foo which, when evaluated,
%   will return a smooth blend between functions provided in handles f1 and f2.
% 
%   The blend is centered around LOC, and most of the smooth blend happens in
%   the region LOC +/- DIST such that F will have value ca. 0.88*f1(x) + 0.12*f2
%   at LOC - DIST and 0.12*f1(x) + 0.88*f2(x) at LOC + DIST. Default LOC = 0.
% 
%   For a pure piecewise function, set the blend distance DIST to zero. F will
%   evaluate to NaN at x = LOC. To avoid evaluating to NaN at x = LOC, use eps
%   instead of zero. Default DIST = 0.
% 
%   f1 and/or f2 may be numeric constants instead of function handles.
% 
%   Example 1: Make a single function handle for a piecewise function (no
%   smoothing).
%     foo = blend(@sin, 1, pi/2);
%     ezplot(foo, [0 pi])
% 
%   Example 2: Make a smooth function from a discontinuous piecewise function
%   and use vectorization to illustrate the effect of different blending
%   distances.
%     f1 = @(x) -x; 
%     f2 = @(x)  x;
%     [dist, x] = meshgrid([0 .01 .1 .5 1], -2:.001:4);
%     foo = blend(f1, f2, 1, dist);
%     plot(x, foo(x))
% 
%   See also tanh.
% 
%   F = BLEND(f1,f2,LOC,DIST)

% Copyright 2016 Sky Sartorius
% Contact: www.mathworks.com/matlabcentral/fileexchange/authors/101715

if nargin < 4
    distance = 0;
end

if nargin < 3 || isempty(location)
    location = 0;
end

validateattributes(location, {'numeric','DimVar'}, {},...
    'blend', 'blending center location', 3);

validateattributes(distance, {'numeric','DimVar'}, {'nonnegative'},...
    'blend', 'blending distance', 4);

if isnumeric(f1)
    f1 = @(x) f1;
end
if isnumeric(f2)
    f2 = @(x) f2;
end

blf = @(x) tanh((x-location)./distance)/2;

foo = @(x) (1/2 - blf(x)).*f1(x) + (1/2 + blf(x)).*f2(x);

% Revision History
%{
2016-04-12 Created.
2016-04-20 
-Added numeric f1 and f2. Creation of function handles may be slower,
but the code is much cleaner than having several versions of the foo = line.
-Simplified and better vectorized based on 1/0 evaluating to inf.
-Input checking.
%}