function valid = timewin_isvalid(T)
% L = TIMEWIN_ISVALID(T) check that the time window is valid/allowable
%       T should be an ordered, disjoint set of N closed intervals,
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

goodsz = size(T,2) == 2;
long = T';
long = long(:);
% comparison should be strict (i.e. not allowing equal values)
goodsort = issorted(long,'strictascend');

valid = goodsz & goodsort;
