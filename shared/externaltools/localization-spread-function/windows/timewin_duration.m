function l = timewin_duration(T)
% L = TIMEWIN_DURATION(T) find the duration of a temporal window
%       T is an ordered, disjoint set of N closed intervals,
%       represented as an Nx2 matrix, so that T(i,1) and T(i,2)
%       are the start and end times of the ith interval, respectively

% Copyright (C) 2021 Thomas Shaw, and Sarah Veatch
% This file is part of SMLM SPACETIME RESOLUTION
% SMLM SPACETIME RESOLUTION is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% SMLM SPACETIME RESOLUTION is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% You should have received a copy of the GNU General Public License
% along with SMLM SPACETIME RESOLUTION.  If not, see <https://www.gnu.org/licenses/>

% check that T is a valid time window
if ~timewin_isvalid(T)
    error('timewin_duration: an invalid time window T was provided');
end

% differences along second dimension, so (end - start) of each interval
% then just add 'em up.
l = sum(diff(T,1,2));
