function [g,errs,time_edge_cor,N,Norm] = spacetime_xcor(x1,y1,t1,x2,y2,t2,spacewin,timewin,...
                                     r,tau,NmVal)
% SPACETIME_XCOR compute spatiotemporal cross correlation function of point data
% [G,ERR,TIME_EDGE_COR,N,NORM] = SPACETIME_XCOR(X1,Y1,T1,X2,Y2,T2,SPACEWIN,TIMEWIN,R,TAU)
%       space-time cross-correlation function of the points X1,Y1,T1, with the
%       points X2,Y2,T2 at R and TAU separations in time and space
%       respectively. SPACEWIN and TIMEWIN specify the spatial window (ROI) and
%       temporal extent of the data, respectively.  Note that R and TAU must be
%       equally spaced (by DR and DTAU respectively) for computational reasons,
%       and that in particular R(1) = DR/2, so that the lower edge of the first
%       r bin is 0.
% [_] = SPACETIME_XCOR(X1,Y1,T1,X2,Y2,T2,SPACEWIN,TIMEWIN,'REdges',R_Edges,'TauEdges',Tau_edges)
%       instead of specifying bin centers R and TAU, the user may specify bin
%       edges.  Bins must still satisfy the same conditions as above (though
%       that may be relaxed in a later release).
% [_] = SPACETIME_XCOR(_,'How', 'Actual')
% [_] = SPACETIME_XCOR(_,'How', 'Uniform') Optional argument 'How' specifies
%       how to do the temporal edge correction, i.e. whether it should assume
%       the observed density or a uniform density in time, respectively.
%       'Actual' is the default.

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

arguments
    x1          (1,:)   double
    y1          (1,:)   double
    t1          (1,:)   double
    x2          (1,:)   double
    y2          (1,:)   double
    t2          (1,:)   double
    spacewin    (1,1)   struct {spacewin_isvalid}
    timewin     (:,2)   double {timewin_isvalid}
    r           (1,:)   double = []
    tau         (1,:)   double = []
    NmVal.REdges    (1,:)   double = []
    NmVal.TauEdges  (1,:)   double = []
    NmVal.How       (1,:)   string = 'actual'
end

    % check that x,y, and t are same size
    if ~isequal(size(x1),size(y1)) || ~isequal(size(x1),size(t1)) || ...
        ~isequal(size(x2),size(y2)) || ~isequal(size(x2),size(t2))
        error('spacetime_xcor: x,y and t must all be the same size')
    end

    % check that the points are in the spatial window
    ind = spacewin_isinside(x1,y1,spacewin) & timewin_isinside(t1,timewin);
    if sum(ind) < numel(ind)
        fprintf(['spacetime_xcor: removing %d points (%.0f %%) that were ',...
            'outside of the ROI\n'], numel(ind) - sum(ind), (1 - sum(ind)/numel(ind))*100);
    end
    x1 = x1(ind); y1 = y1(ind); t1 = t1(ind);
    ind = spacewin_isinside(x2,y2,spacewin) & timewin_isinside(t2,timewin);
    if sum(ind) < numel(ind)
        fprintf(['spacetime_xcor: removing %d points (%.0f %%) that were ',...
            'outside of the ROI\n'], numel(ind) - sum(ind), (1 - sum(ind)/numel(ind))*100);
    end
    x2 = x2(ind); y2 = y2(ind); t2 = t2(ind);

    % Check that r and tau satisfy restrictions
    rbinedges = NmVal.REdges;
    if isempty(r) && isempty(rbinedges)
        error('spacetime_xcor: no r values were specified')
    elseif isempty(r)
        r = rbinedges(2:end) - diff(rbinedges)/2;
    elseif isempty(rbinedges)
        dr = r(2) - r(1);
        rbinedges = [r - dr/2, r(end) + dr/2];
    elseif ~isequal(r, rbinedges(2:end) - diff(rbinedges)/2)
        error('spacetime_xcor: incompatible r and REdges')
    end

    taubinedges = NmVal.TauEdges;
    if isempty(tau) && isempty(taubinedges)
        error('spacetime_xcor: no tau values were specified')
    elseif isempty(tau)
        tau = taubinedges(2:end) - diff(taubinedges)/2;
    elseif isempty(taubinedges)
        dtau = tau(2) - tau(1);
        taubinedges = [tau - dtau/2, tau(end) + dtau/2];
    elseif ~isequal(tau, taubinedges(2:end) - diff(taubinedges)/2)
        error('spacetime_xcor: incompatible tau and TauEdges')
    end
        
    difftau = diff(tau);
    if (max(difftau) - min(difftau))/min(difftau) > 1e-13
        error('spacetime_xcor: requested tau values must be equally spaced. Support for unequally spaced tau may be added in a future release.')
    end
    Dtau = difftau(1);
    taubinedges = min(tau)-Dtau/2 : Dtau : max(tau)+Dtau/2;

    diffr = diff(r);
    Dr = diffr(1);
    if (max(diffr) - min(diffr))/min(diffr) > 1e-13 
        fprintf('Smallest r bin width: %f, largest %f\n', min(diffr), max(diffr));
        error('spactime_xcor: requested r values must be equally spaced. Support for unequally spaced r may be added in a future release.')
    elseif abs((r(1) - Dr/2)/Dr) > 1e-13
        fprintf('First r bin: center: %f, lower edge: %f\n', r(1), rbinedges(1));
        error('spacetime_xcor: , and smallest r bin must start at 0')
    end
    % clean them up
    rbinedges = (0:numel(r))*Dr;
    r = rbinedges(2:end) - Dr/2;
    
    taumin = min(taubinedges);
    taumax = max(taubinedges);
    rmax = max(rbinedges);
    
    % N is just the histogram of pairs, in r and tau bins
    N = crosspairs_st_binned(x1,y1,t1, x2,y2,t2, rmax, numel(r), taumin, taumax, numel(tau));
    
    % basic normalization for area and density (no edge corrections)
    area_per_rbin = 2*pi*r'*Dr;
    time_per_tbin = Dtau;
    area = spacewin_area(spacewin);
    duration_excluding_gaps = timewin_duration(timewin);
    
    density1 = numel(x1)/area/duration_excluding_gaps;
    density2 = numel(x2)/area/duration_excluding_gaps;
    
    basic_normalization = duration_excluding_gaps*area*density1*density2*area_per_rbin*time_per_tbin;
    
    % edge corrections: space first, then time
    edge_cor = spatial_edge_correction(spacewin, r);

    how = lower(NmVal.How);
    if strcmp(how, 'uniform')
        time_edge_cor = time_edge_correction_unif(taubinedges, timevec);
    elseif strcmp(how, 'actual')
        time_edge_cor = time_edge_correction_density_cross(t1,t2,taubinedges,timewin);
    else
        error('invalid time edge correction method supplied')
    end
    
    % Norm is the normalization that turns N into an edge-corrected correlation function
    Norm = basic_normalization.*time_edge_cor.*edge_cor;
    g = N./Norm;
    
    errs = sqrt(N)./Norm;

end
