function inds = timewin_isinside(t, timewin)
% TIMEWIN_ISINSIDE check which points are in a given temporal window
% INDS = TIMEWIN_ISINSIDE(T,TIMEWIN)    Return a logical vector the same size as T
%           such that INDS(i) is true iff T(i) is inside TIMEWIN

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

if ~timewin_isvalid(timewin)
    error('timewin_isinside: invalid temporal window provided');
end

% t must be sorted to enable a fast algorithm
[sortedt,sorting] = sort(t);

% initialize vars for keeping track of stuff
jfirst = find(timewin(:,1) <= sortedt(1),1); % index into twflat
jlast = find(timewin(:,2) > sortedt(end),1,'last');
if isempty(jfirst)
    jfirst = 1;
end
if isempty(jlast)
    jlast = size(timewin,1);
end

inds = false(size(sortedt));
for j = jfirst:jlast
    ifirst = find(sortedt >= timewin(j,1),1);
    ilast = find(sortedt < timewin(j,2),1,'last');

    inds(ifirst:ilast) = true;
end

reorderedinds = sorting(inds);
inds = false(size(t));
inds(reorderedinds) = true;
end
