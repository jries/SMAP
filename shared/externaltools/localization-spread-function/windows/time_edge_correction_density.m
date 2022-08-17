function [taufactor,norm] = time_edge_correction_density(t,tau_edges,timewin)
% TIME_EDGE_CORRECTION_DENSITY temporal edge correction for correlation
% function

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

ntout = numel(tau_edges) -1;
% check validity of timewin argument
if ~timewin_isvalid(timewin)
    error('time_edge_correction: invalid time window provided');
end

timevec = sort(unique(t));
%fprintf('There are %f times as many points as unique times\n', numel(t)/numel(timevec));

Nperframe = arrayfun(@(tt) sum(t == tt), timevec);

tmax = numel(timevec);

timediffs = zeros(tmax*(tmax-1)/2, 1);
weights = zeros(tmax*(tmax-1)/2, 1);
count = 1;
for i = 1:tmax-1
    for j = i:tmax
        timediffs(count) = timevec(j) - timevec(i);
        weights(count) = Nperframe(i)*Nperframe(j);
        count = count + 1;
    end
end

dtau = diff(tau_edges);

[~, ~, bin] = histcounts(timediffs, tau_edges);
inds = bin>0;
exptauperbin = accumarray(bin(inds), weights(inds)',[ntout,1])./dtau(:);

norm = numel(t)^2 / timewin_duration(timewin);
taufactor = exptauperbin'/norm;
