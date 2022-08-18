function [out, params] = spacetime_resolution(varargin)
% SPACETIME_RESOLUTION compute spatiotemporal point-spread function from point data
% OUT = SPACETIME_RESOLUTION(X,Y,T,SPACEWIN,TIMEWIN)    Compute the spacetime
%               resolution for the points given by X,Y,T (vectors with a value for each point)
%               on the spatial window (ROI) SPACEWIN, (as validated by spacewin_isvalid()), and the
%               temporal window TIMEWIN (as validated by timewin_isvalid())
% OUT = SPACETIME_RESOLUTION(DATA)      As above, but with the arguments given as fields of a struct,
%               with fieldnames 'x','y','t','spacewin', and 'timewin'.
% OUT = SPACETIME_RESOLUTION(_,NAME, VALUE)     Optional arguments may be given as name-value pairs:
%               'R' (default: 2.5:5:250): vector of spatial separation values to use for the computation.
%                               It may be appropriate to change this if your data is given in units other
%                               than nm.
%               'TauBinMethod' (default: 'logspace'): string specifying the method to use for automatic
%                               determination of tau ranges to calculate resolution for.
%               'NTauBin' (default: 10): Number of tau bins to calculate resolution for.
%               'TauEdges' (no default): Explicit tau bin edges. Overrides TauBinMethod and NTauBin.
%               'SigmaStartPt' (default: 10): Starting point for resolution in nonlinear fitting step.
%                               It may be appropriate to change this if your data is given in units other
%                               than nm.
%               'Bootstrap' (default: false): Generate a bootstrapped confidence interval by resampling the
%                               data. This is off by default because it takes a considerable amount of time.

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

%% Deal with input arguments and parameters
args = varargin;
[data, params] = handle_args(args); % defined below

r = params.r;

x = data.x;
y = data.y;
t = data.t;
spacewin = data.spacewin;
timewin = data.timewin;
% only use points that are inside the spatial window
inds = spacewin_isinside(x,y,spacewin) & timewin_isinside(t,timewin);
if ~all(inds)
    fprintf('spatial and temporal windows contain %d of %d points (%.0f %%)\n',...
        sum(inds), numel(inds), sum(inds)./numel(inds) * 100);
    x = x(inds);
    y = y(inds);
    t = t(inds);
end

timevec = sort(unique(t));
% TODO: optionally specify frame times... somehow?
frame_time = timevec(2) - timevec(1);
taumax = timevec(end) - timevec(1);
tau = frame_time:frame_time:taumax;

[g_alltau, ~, time_edge_cor, counts_all, norm_all] = spacetime_acor(...
                    x,y,t,spacewin,timewin,r,tau);

%% Choose larger tau bins if applicable

switch params.TauBinMethod
    case 'twowidth'
        [~, taumaxind] = guess_taumax_from_g(g_alltau, r, tau, 25, .5);
        bins = guess_tau_edges_from_time_edge_cor(time_edge_cor(1:taumaxind), params.NTauBin);
        bins(end+1) = numel(tau);
        taubinedges = tau(bins);
        taubincenters = taubinedges(2:end) - diff(taubinedges)/2;
    case 'logspace'
        [~, taumaxind] = guess_taumax_from_g(g_alltau, r, tau, 25, .5);
        bins = guess_tau_logspace(taumaxind, params.NTauBin);
        bins(end+1) = numel(tau);
        taubinedges = tau(bins);
        taubincenters = taubinedges(2:end) - diff(taubinedges)/2;
    case 'provided'
        taubinedges = params.TauEdges;
        for i =1:numel(taubinedges)
            bins(i) = find(tau >= taubinedges(i),1);
        end
        taubincenters = taubinedges(2:end) - diff(taubinedges)/2;
end

%% average g(r,tau) into the larger tau bins
% this takes a weighted average over the tau bins (might not matter)
norm_g = zeros(numel(r),numel(taubincenters));
counts_g = zeros(size(norm_g));
for i=1:numel(bins)-1
    inds = bins(i):bins(i+1)-1;
    gg = g_alltau(:, inds);
    
    norm_g(:,i) = sum(norm_all(:,inds),2);
    counts_g(:,i) = sum(counts_all(:,inds),2);
    dg(:, i) = std(gg, [], 2);
end

g = counts_g./norm_g;

%% compute a Coefficient of Variation based on an assumption of poisson.
rinds = 1:find(r > 2*params.SigmaStartpoint,1);
Ns = sum(counts_g(rinds,:),1);
Norms = sum(norm_g(rinds,:),1);
COV = sqrt(Ns(1:end-1)./Norms(1:end-1).^2 + Ns(end)./Norms(end).^2)./(Ns(1:end-1)./Norms(1:end-1) - Ns(end)./Norms(end));

%% subtract the long tau bin
g_diff = g(:, 1:end-1)-g(:, end);

% normalize to simplify fitting
g_diff_norm = g_diff./g_diff(1, :);

%% fit the subtracted curves
%sigma startpoint is an optional parameter
startpoint = [1 params.SigmaStartpoint];
fitpts = r<100;
for i=1:size(g_diff_norm, 2)
    F = fit(r(fitpts)', g_diff_norm(fitpts, i), 'A*exp(-x.^2/4/s^2)', 'startpoint', startpoint, 'lower', [0 5]);
    s(i) = F.s;
    ci = diff(confint(F, .68));
    
    ds(i) = ci(2)/2;
end

minds = .5;
%weights = 1./max(minds, ds).^2;
weights = sum(counts_g(:,1:end),1); % all r in given tau bin
S = sum(weights.*[s s(end)])/sum(weights);
dS = sqrt(sum([ds, ds(end)].^2.*weights.^2))/sum(weights);

%% bootstrap CIs
conf = nan(size(s));
if params.Bootstrap
    fprintf('spacetime_resolution: starting bootstrap for confidence intervals\n');
    tic;
    npts = numel(x);
    npts_sample = round(npts/4);
    s_resamp = zeros(8,numel(s));
    for i=1:8
        inds = datasample(1:npts, npts_sample, 'replace', false);
        [~,~,~,countsi, normsi] = spacetime_acor(x(inds), y(inds),t(inds),spacewin,timewin,r,tau);
        [gi] = reaverage_grtau(countsi,normsi,bins);
        dg = gi(:,1:end-1) - gi(:,end);
        dg_norm = dg ./ dg(1,:);
        for j=1:size(dg_norm, 2)
            F = fit(r(fitpts)', dg_norm(fitpts, j), 'A*exp(-x.^2/4/s^2)', 'startpoint', startpoint, 'lower', [0 5]);
            s_resamp(i,j) = F.s;
        end
        fprintf('%d ',i);
    end
    fprintf('\n');
    conf = std(s_resamp,[],1)/2;
    dS = sqrt(sum([conf, conf(end)].^2.*weights.^2))/sum(weights);

    toc
end

%% Assemble outputs
out = struct('S', S, 'dS', dS, 's',s, 'ds', ds, 'Dg', g_diff, 'nDg', g_diff_norm,...
    'taubincenters', taubincenters, 'cWA', g, 'dcWA', dg, 'taubinedges', taubinedges,...
    'c', g_alltau, 'COV', COV, 'confint', conf,'r', r, 'tau', tau);
end

function [g,counts, norm] = reaverage_grtau(counts_all, norm_all, taubin_inds)
    counts = zeros(size(counts_all,1),numel(taubin_inds)-1);
    norm = zeros(size(counts));
    for i = 1:numel(taubin_inds)-1
        inds = taubin_inds(i):taubin_inds(i+1)-1;
        norm(:,i) = sum(norm_all(:,inds),2);
        counts(:,i) = sum(counts_all(:,inds),2);
    end
    
    g = counts./norm;
end


function [taumax, taumaxind] = guess_taumax_from_g(g, r, tau, rmax, threshold)
inds = r<=rmax;
dr = r(2)-r(1);
gvstau_shortr = 2*pi*dr*sum(r(inds)'.*g(inds, :))/pi/rmax^2;
taumaxind = find(gvstau_shortr/mean(gvstau_shortr(round(end/4*3):end))<(threshold+1), 1);
if isempty(taumaxind) || taumaxind > numel(tau)/4*3
    taumaxind = round(numel(tau)/4*3); 
end
taumax = tau(taumaxind);
end

function tau_bins = guess_tau_logspace(ntau,nbin)
tau_bins = unique(round(logspace(log10(1), log10(ntau), nbin)));
%tau_bins(end+1) = ntau;
end


function tau_bins = guess_tau_edges_from_time_edge_cor(time_edge_cor,nbin)
totpairs = cumsum(time_edge_cor);
nlargebins = floor(nbin/2);
nsmallbins = nbin - nlargebins;
binratio = .2;

N = max(totpairs);

dd = N/(1+nlargebins+nsmallbins*binratio);
pairedges = [(0:nsmallbins-1)*dd*binratio nsmallbins*dd*binratio:dd:N];

tau_bins =ones(size(pairedges));

for i=2:numel(pairedges)-1
    tau_bins(i) = find(totpairs>pairedges(i), 1);
end
tau_bins(end) = numel(totpairs);
end

function [tau_bins, tau_edges, tau_centers] = guess_tau_equal_from_cumsum(metric, tau, nbin)
cs = cumsum(metric);
binsize = cs(end)/nbin;


for i=1:nbin
tau_bins(i) = find(cs > (i-1)*binsize,1);
end
if tau_bins(end) ~= tau(end)
    tau_bins = [tau_bins,tau(end)];
end
tau_edges=tau_bins;
tau_centers = tau_edges(2:end) - diff(tau_edges)/2;
end

function [data, params] = handle_args(args)
% data has the point process data itself:
% x,y,t,spacewin,timewin
% params has extra parameters

if numel(args) == 0
    error('spacetime_resolution: no data supplied')
end

% handle args passed as struct
argnames = {'x', 'y', 't', 'spacewin', 'timewin'};
if isstruct(args{1})
    if all(ismember(argnames, fieldnames(args{1})))
        dd = args{1};
        args = args(2:end);
        for i=1:numel(argnames)
            f = argnames{i};
            data.(f) = dd.(f);
        end
    else
        error('spacetime_resolution: DATA must contain fields x, y, t, spacewin, and timewin');
    end
else
    % no data struct provided
    if numel(args) < 5
        error('spacetime_resolution: not enough data provided')
    else
        % assume x,y,t,spacewin,timewin are the next five
        [data.x,data.y,data.t,data.spacewin,data.timewin] = deal(args{1:5});
        args = args(6:end);
    end
end
% have data at this point
sz = size(data.x);
if ~isequal(sz, size(data.y)) || ~isequal(sz, size(data.t))
    error('spacetime_resolution: sizes of x,y, and t must match');
end
if ~timewin_isvalid(data.timewin)
    error('spacetime_resolution: timewin is not valid')
end
if ~spacewin_isvalid(data.spacewin)
    error('spacetime_resolution: spacewin is not valid')
end

%% handle remaining parameters
% Set defaults
rdefault = 2.5:5:250;
taubinmethods = {'twowidth', 'logspace', 'provided'};

p = inputParser;
p.FunctionName = 'spacetime_resolution';
p.StructExpand = true;
p.KeepUnmatched = true;

p.addParameter('r', rdefault, @isnumeric);
p.addParameter('TauBinMethod', 'logspace', @(x) any(validatestring(x,taubinmethods)));
p.addParameter('NTauBin', 10, @(x) isnumeric(x) && isscalar(x))
p.addParameter('TauEdges', [], @(x) isnumeric(x) && issorted(x))
p.addParameter('SigmaStartpoint', 10, @(x) isnumeric(x) && isscalar(x))
p.addParameter('Bootstrap', false)
%TODO: add more

p.parse(args{:});
params = p.Results;

%warn if there are extras
if numel(fieldnames(p.Unmatched)) > 0
    warning('spacetime_resolution: some supplied parameters were not recognized. See below:');
    disp(p.Unmatched);
end

% normalize conflicting options
if ~isempty(params.TauEdges)
    if ~any(strcmp('TauBinMethod', p.UsingDefaults)) || ~any(strcmp('NTauBin',p.UsingDefaults))
        warning('spacetime_resolution: setting TauEdges overrides NTauBin and TauBinMethod')
    end
    params.TauBinMethod = 'provided';
    params.NTauBin = numel(params.TauEdges) - 1;
end
end
