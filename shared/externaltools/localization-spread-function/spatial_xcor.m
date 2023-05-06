function [g,errs,N,Norm] = spatial_xcor(x1,y1,x2,y2,spacewin,r,NmVal)
% SPATIAL_XCOR compute spatial cross correlation of two point data sets
% [G,ERR,N,NORM] = SPATIAL_XCOR(X1,Y1,X2,Y2,T2,SPACEWIN,R)
%       spatial cross-correlation function of the points X1,Y1, with the points
%       X2,Y2 at R separations in space. SPACEWIN specifies the spatial window
%       (ROI) for the data. Note that R must be equally spaced (by DR) for
%       computational reasons, and that in particular R(1) = DR/2, so that the
%       lower edge of the first r bin is 0.
%       SPATIAL_XCOR(X1,Y1,X2,Y2,SPACEWIN,'REdges',R_Edges)
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
    x1          (1,:)   double
    y1          (1,:)   double
    x2          (1,:)   double
    y2          (1,:)   double
    spacewin    (1,1)   struct {spacewin_isvalid}
    r           (1,:)   double = []
    NmVal.REdges    (1,:)   double = []
end

    % check that x and y are same size
    if ~isequal(size(x1),size(y1)) || ~isequal(size(x2),size(y2))
        error('spatial_xcor: x and y must be the same size')
    end

    % check that the points are in the spatial window
    ind = spacewin_isinside(x1,y1,spacewin);
    if sum(ind) < numel(ind)
        fprintf(['spatial_xcor: removing %d points (%.0f %%) that were ',...
            'outside of the ROI\n'], numel(ind) - sum(ind), (1 - sum(ind)/numel(ind))*100);
    end
    x1 = x1(ind); y1 = y1(ind);
    ind = spacewin_isinside(x2,y2,spacewin);
    if sum(ind) < numel(ind)
        fprintf(['spatial_xcor: removing %d points (%.0f %%) that were ',...
            'outside of the ROI\n'], numel(ind) - sum(ind), (1 - sum(ind)/numel(ind))*100);
    end
    x2 = x2(ind); y2 = y2(ind);

    % Check that r and tau satisfy restrictions
    rbinedges = NmVal.REdges;
    if isempty(r) && isempty(rbinedges)
        error('spatial_xcor: no r values were specified')
    elseif isempty(r)
        r = rbinedges(2:end) - diff(rbinedges)/2;
    elseif isempty(rbinedges)
        dr = r(2) - r(1);
        rbinedges = [r - dr/2, r(end) + dr/2];
    elseif ~isequal(r, rbinedges(2:end) - diff(rbinedges)/2)
        error('spatial_xcor: incompatible r and REdges')
    end

    diffr = diff(r);
    Dr = diffr(1);
    if (max(diffr) - min(diffr))/min(diffr) > 1e-13 
        fprintf('Smallest r bin width: %f, largest %f\n', min(diffr), max(diffr));
        error('spatial_xcor: requested r values must be equally spaced. Support for unequally spaced r may be added in a future release.')
    elseif abs((r(1) - Dr/2)/Dr) > 1e-13
        fprintf('First r bin: center: %f, lower edge: %f\n', r(1), rbinedges(1));
        error('spatial_xcor: , and smallest r bin must start at 0')
    end
    % clean them up
    rbinedges = (0:numel(r))*Dr;
    r = rbinedges(2:end) - Dr/2;
    
    rmax = max(rbinedges);
    
    % N is just the histogram of pairs, in r and tau bins
    N = crosspairs_binned(x1,y1, x2,y2, rmax, numel(r));
    
    % basic normalization for area and density (no edge corrections)
    area_per_rbin = 2*pi*r'*Dr;
    area = spacewin_area(spacewin);
    
    density1 = numel(x1)/area;
    density2 = numel(x2)/area;
    
    basic_normalization = area*density1*density2*area_per_rbin;
    
    % edge corrections
    edge_cor = spatial_edge_correction(spacewin, r);

    % Norm is the normalization that turns N into an edge-corrected correlation function
    Norm = basic_normalization.*edge_cor;
    g = N./Norm;
    
    errs = sqrt(N)./Norm;

end
