function [g,errs,time_edge_cor,N,Norm] = spacetime_acor_xy(x,y,t,spacewin,timewin,hxy,tau,NmVal)
% SPACETIME_ACOR_XY compute spatiotemporal correlation function for point data
% [G,ERR,TIME_EDGE_COR,N,NORM] = SPACETIME_ACOR_XY(X,Y,T,SPACEWIN,TIMEWIN,R,TAU)
%       space-time autocorrelation function of the points X,Y,T, at a lattice
%       of points specified by HXY in x and y and TAU separations in time and
%       space respectively. SPACEWIN and TIMEWIN specify the spatial window
%       (ROI) and temporal extent of the data, respectively. HXY and TAU should
%       be equally spaced.
% [_] = SPACETIME_ACOR_XY(X,Y,T,SPACEWIN,TIMEWIN,'XYEdges',XY_Edges,'TauEdges',Tau_edges)
%       instead of specifying bin centers HXY and TAU, the user may specify bin edges.
% [_] = SPACETIME_ACOR_XY(_,'How', 'Actual')
% [_] = SPACETIME_ACOR_XY(_,'How', 'Uniform') Optional argument 'How' specifies how to do
%       the temporal edge correction, i.e. whether it should assume the observed density
%       or a uniform density in time, respectively. 'Actual' is the default.

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
    t           (1,:)   double
    spacewin    (1,1)   struct {spacewin_isvalid}
    timewin     (:,2)   double {timewin_isvalid}
    hxy           (1,:)   double = []
    tau         (1,:)   double = []
    NmVal.XYEdges    (1,:)   double = []
    NmVal.TauEdges  (1,:)   double = []
    NmVal.How       (1,:)   string = 'actual'
end

    % check that x,y, and t are same size
    if ~isequal(size(x),size(y)) || ~isequal(size(x),size(t))
        error('spacetime_acor: x,y and t must all be the same size')
    end

    % check that the points are in the spatial window
    ind = spacewin_isinside(x,y,spacewin) & timewin_isinside(t,timewin);
    if sum(ind) < numel(ind)
        fprintf(['spacetime_acor: removing %d points (%.0f %%) that were ',...
            'outside of the ROI or temporal window\n'], numel(ind) - sum(ind), (1 - sum(ind)/numel(ind))*100);
    end
    x = x(ind); y = y(ind); t = t(ind);

    % Check that r and tau satisfy restrictions
    xyedges = NmVal.XYEdges;
    if isempty(hxy) && isempty(xyedges)
        error('spacetime_acor: no hxy values were specified')
    elseif isempty(hxy)
        hxy = xyedges(2:end) - diff(xyedges)/2;
    elseif isempty(xyedges)
        dxy = hxy(2) - hxy(1);
        xyedges = [hxy - dxy/2, hxy(end) + dxy/2];
    elseif ~isequal(hxy, xyedges(2:end) - diff(xyedges)/2)
        error('spacetime_acor: incompatible hxy and XYEdges')
    end

    taubinedges = NmVal.TauEdges;
    if isempty(tau) && isempty(taubinedges)
        error('spacetime_acor: no tau values were specified')
    elseif isempty(tau)
        tau = taubinedges(2:end) - diff(taubinedges)/2;
    elseif isempty(taubinedges)
        dtau = tau(2) - tau(1);
        taubinedges = [tau - dtau/2, tau(end) + dtau/2];
    elseif ~isequal(tau, taubinedges(2:end) - diff(taubinedges)/2)
        error('spacetime_acor: incompatible tau and TauEdges')
    end

    nxy = numel(hxy);
    dxy = diff(xyedges);
    ntau = numel(tau);
    dtau = diff(taubinedges);

    taumin = min(taubinedges);
    taumax = max(taubinedges);
    rmax = max(abs(xyedges))*sqrt(2);
    
    % N is just the histogram of pairs, in r and tau bins
    noutmax = 2e8;
    [dx,dy,dt,err] = closepairs_st(x,y,t, rmax, taumin, taumax,noutmax);
    if err
        error('closepairs_st: number of pairs exceeds noutmax. Try using a smaller dataset or increasing noutmax')
    end
    % quick histcounts3
    ind = zeros(numel(dx),3);
    [~,ind(:,1)] = histc(dy, xyedges); % first dimension is y, following matlab convention
    [~,ind(:,2)] = histc(dx, xyedges);
    [~,ind(:,3)] = histc(dt, taubinedges);
    valid = all(ind > 0, 2); % exclude points that are not in any bins
    N = accumarray(ind(valid,:), 1, [nxy, nxy, ntau]);
    
    % basic normalization for area and density (no edge corrections)
    area_per_xybin = dxy.*dxy';
    time_per_tbin = reshape(dtau,1,1,numel(dtau));
    area = spacewin_area(spacewin);
    duration_excluding_gaps = timewin_duration(timewin);
    
    density = numel(x)/area/duration_excluding_gaps;
    
    basic_normalization = duration_excluding_gaps*area*density^2.*area_per_xybin.*time_per_tbin;
    
    % edge corrections: space first, then time
    [hx, hy] = meshgrid(hxy, hxy);
    edge_cor = spatial_edge_correction(spacewin, hx, hy);

    how = lower(NmVal.How);
    if strcmp(how, 'uniform')
        time_edge_cor = time_edge_correction_unif(taubinedges, timevec);
    elseif strcmp(how, 'actual')
        time_edge_cor = time_edge_correction_density(t,taubinedges,timewin);
    else
        error('invalid time edge correction method supplied')
    end
    time_edge_cor = reshape(time_edge_cor, 1,1, numel(time_edge_cor));
    
    % Norm is the normalization that turns N into an edge-corrected correlation function
    Norm = basic_normalization.*time_edge_cor.*edge_cor;
    g = N./Norm;
    
    errs = sqrt(N)./Norm;

end
