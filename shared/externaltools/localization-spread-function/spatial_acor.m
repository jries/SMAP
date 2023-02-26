function [g,errs,N,Norm] = spatial_acor(x,y,spacewin,r,NmVal)
% SPATIAL_ACOR compute spatial correlation function of point data
% [G,ERR,N,NORM] = SPATIAL_ACOR(X,Y,SPACEWIN,R)
%       space-time autocorrelation function of the points X,Y, at R separations
%       in space. SPACEWIN specifies the spatial window (ROI) of the data.
%       Note that R must be equally spaced (by DR) for computational reasons,
%       and that in particular R(1) = DR/2, so that the lower edge of the first
%       r bin is 0.
%       SPATIAL_ACOR(X,Y,SPACEWIN,'REdges',R_Edges)
%       instead of specifying bin centers R, the user may specify bin edges.
%       Bins must still satisfy the same conditions as above (though that may
%       be relaxed in a later release).

% Copyright (C) 2022 Thomas Shaw, and Sarah Veatch
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

arguments
    x           (1,:)   double
    y           (1,:)   double
    spacewin    (1,1)   struct {spacewin_isvalid}
    r           (1,:)   double = []
    NmVal.REdges    (1,:)   double = []
end

    % check that x,y are same size
    if ~isequal(size(x),size(y))
        error('spatial_acor: x and y must be the same size')
    end

    % check that the points are in the spatial window
    ind = spacewin_isinside(x,y,spacewin);
    if sum(ind) < numel(ind)
        fprintf(['spatial_acor: removing %d points (%.0f %%) that were ',...
            'outside of the ROI\n'], numel(ind) - sum(ind), (1 - sum(ind)/numel(ind))*100);
    end
    x = x(ind); y = y(ind);

    % Check that r satisfy restrictions
    rbinedges = NmVal.REdges;
    if isempty(r) && isempty(rbinedges)
        error('spatial_acor: no r values were specified')
    elseif isempty(r)
        r = rbinedges(2:end) - diff(rbinedges)/2;
    elseif isempty(rbinedges)
        dr = r(2) - r(1);
        rbinedges = [r - dr/2, r(end) + dr/2];
    elseif ~isequal(r, rbinedges(2:end) - diff(rbinedges)/2)
        error('spatial_acor: incompatible r and REdges')
    end

    diffr = diff(r);
    Dr = diffr(1);
    if (max(diffr) - min(diffr))/min(diffr) > 1e-13 
        fprintf('Smallest r bin width: %f, largest %f\n', min(diffr), max(diffr));
        error('spatial_acor: requested r values must be equally spaced. Support for unequally spaced r may be added in a future release.')
    elseif abs((r(1) - Dr/2)/Dr) > 1e-13
        fprintf('First r bin: center: %f, lower edge: %f\n', r(1), rbinedges(1));
        error('spatial_acor: , and smallest r bin must start at 0')
    end
    % clean them up
    rbinedges = (0:numel(r))*Dr;
    r = rbinedges(2:end) - Dr/2;
    
    rmax = max(rbinedges);
    
    % N is just the histogram of pairs, in r and tau bins
    N = closepairs_binned(x,y, rmax, numel(r));
    
    % basic normalization for area and density (no edge corrections)
    area_per_rbin = 2*pi*r'*Dr;
    area = spacewin_area(spacewin);
    
    density = numel(x)/area;
    
    basic_normalization = area*density^2*area_per_rbin;
    
    % edge corrections
    edge_cor = spatial_edge_correction(spacewin, r);

    % Norm is the normalization that turns N into an edge-corrected correlation function
    Norm = basic_normalization.*edge_cor;
    g = N./Norm;
    
    errs = sqrt(N)./Norm;

end
