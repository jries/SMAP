function taufactor = time_edge_correction_unif(tau_edges, timewin)
% TIME_EDGE_CORRECTION_UNIF temporal edge correction assuming uniform rate

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

tau = tau_edges(2:end) - diff(tau_edges)/2;

if ~timewin_isvalid(timewin)
    error('time_edge_correction: invalid time window provided');
end

T = timewin_duration(timewin);
overlap = arrayfun(@(t) timewin_overlap(timewin, timewin + t), tau);

taufactor = overlap/T;
end