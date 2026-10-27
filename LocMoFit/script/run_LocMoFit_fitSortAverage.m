function out = run_LocMoFit_fitSortAverage(smlFile, varargin)
% RUN_LOCMOFIT_FITSORTAVERAGE  Fit, sort, align and average the sites of one _sml.mat.
%
%   out = run_LocMoFit_fitSortAverage(smlFile)
%   out = run_LocMoFit_fitSortAverage(smlFile, 'fit', true, 'nBins', 10, ...)
%
% Runs the ROI-manager plugin LocMoFit_fitSortAverage on a SMAP file whose sites
% are marked and rotated, without clicking: it starts SMAP, loads the file,
% restores the file's render layers (LocMoFit fits on them), fits the sites if
% needed, sorts them, averages them in bins aligned on the fitted base, and
% writes next to the input (or to 'outDir'):
%   <name>_sml.mat   the bin averages only, as a SMAP file (load it in SMAP)
%   <name>.png       the averages rendered with the file's render layers, and
%                    <name>_layer<k>.png for each layer alone
%   <name>.csv       the bins: sites, number of sites, median sort value
% The input file is only read. An existing output is never overwritten.
%
% Options (name-value):
%   'settings'   LocMoFit settings file (*_LocMoFit.mat) loaded into the
%                evaluator. Default: the Myo5/Sla2 settings that come with the
%                plugin (plugins/+ROIManager/+Analyze/LocMoFit_settings). Only
%                its models are used when 'fit' is false.
%   'evaluator'  Name of the LocMoFit evaluation in the sites. Default 'LocMoFitGUI'.
%   'fit'        true: fit all sites; false: use the saved fits; 'auto'
%                (default): fit only if no site with use checked has a fit yet.
%                Sites without a fit are otherwise left out (and counted).
%                Fitting takes about 16 s per site.
%   'sortPar'    Fit parameter to sort by. Default 'm2.lPar.y' (base-cap distance).
%   'direction'  'descend' (default: flat sites in bin 1, left) or 'ascend'.
%   'selection'  'use' (default: sites with use checked) or 'all'.
%   'excludeSites' Site IDs to leave out (their use is unchecked in memory only;
%                the input file is not changed). Default none.
%   'binBy'      'sites' (default): nBins bins with equal numbers of sites.
%                'distance': one bin per range of the sort value between
%                consecutive binEdges, the same ranges for every data set.
%   'nBins'      binBy 'sites': number of bins, from the first to the last
%                usable site. Default 10.
%   'binEdges'   binBy 'distance': bin edges (nm), all decreasing or all
%                increasing; bins appear left to right in this order, each with
%                its lower edge included; sites outside all bins are left out.
%                Default [Inf 15 0 -10 -20 -35 -50 -70 -95 -125 -Inf].
%   'crop'       'site ROI' (default): each site contributes the localizations of
%                its site ROI box, as in AverageSites_window, aligned on the fit;
%                the averages fade out softly at the edge. 'aligned square':
%                every site cut to the same square after alignment (sharp edge).
%   'boxSize'    'aligned square' only: side (nm) of the square. Default: 1.5 x
%                the site ROI size.
%   'removeHotSpots' true (default): leave out places where one site emits again
%                and again (dye stuck on, bead, dirt); see the plugin.
%   'hotFactor'  Threshold factor for hot spots. Default 3.
%   'layout'     'fixed' (default): the fitted base of every site at (b x
%                binSpacing, binSpacing) for bin b, identical in every data set,
%                so averages of different data sets overlay exactly. 'compact':
%                centred on the ROI centres, with binGap between the averages.
%   'binSpacing' Fixed layout: distance (nm) between bin centres. Default 550.
%                Use the same value for all data sets you want to combine.
%   'binGap'     Compact layout: empty space (nm) between neighbouring averages, edge to edge
%                (edge: 95 % of the localizations, about where they visibly
%                end). Default 250.
%   'outDir'     Output folder. Default: the input's folder.
%   'name'       Output name. Default: <input name>_aligned_avg_<number of bins>bins
%                (with '_bydistance' for binBy 'distance').
%
% Example:
%   run_LocMoFit_fitSortAverage('/data/Myo5_marked_sml.mat', 'fit', true)

ip = inputParser;
smapRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));  % SMAP/LocMoFit/script -> SMAP
ip.addParameter('settings', fullfile(smapRoot, 'plugins', '+ROIManager', '+Analyze', 'LocMoFit_settings', 'Myo5_Sla2_settings_ver2_LocMoFit.mat'));
ip.addParameter('evaluator', 'LocMoFitGUI');
ip.addParameter('fit', 'auto');
ip.addParameter('sortPar', 'm2.lPar.y');
ip.addParameter('direction', 'descend');
ip.addParameter('selection', 'use');
ip.addParameter('excludeSites', []);
ip.addParameter('nBins', 10);
ip.addParameter('binBy', 'sites');
ip.addParameter('binEdges', [Inf 15 0 -10 -20 -35 -50 -70 -95 -125 -Inf]);
ip.addParameter('crop', 'site ROI');
ip.addParameter('boxSize', []);
ip.addParameter('binGap', 250);
ip.addParameter('layout', 'fixed');
ip.addParameter('binSpacing', 550);
ip.addParameter('removeHotSpots', true);
ip.addParameter('hotFactor', 3);
ip.addParameter('outDir', '');
ip.addParameter('name', '');
ip.parse(varargin{:});
o = ip.Results;

smlFile = char(smlFile);
[inDir, stem] = fileparts(smlFile);
if isempty(o.outDir), o.outDir = inDir; end
byDistance = strcmp(o.binBy, 'distance');
nBinsOut = o.nBins; if byDistance, nBinsOut = numel(o.binEdges) - 1; end
if isempty(o.name)
    o.name = sprintf('%s_aligned_avg_%dbins', regexprep(stem, '_sml$', ''), nBinsOut);
    if byDistance, o.name = [o.name '_bydistance']; end
end
outSml = fullfile(o.outDir, [o.name '_sml.mat']);
outPng = fullfile(o.outDir, [o.name '.png']);
outCsv = fullfile(o.outDir, [o.name '.csv']);
for f = {outSml, outPng, outCsv}
    if isfile(f{1}), error('%s exists; choose another ''name'' or ''outDir''.', f{1}); end
end
if ~isfile(o.settings), error('LocMoFit settings not found: %s', o.settings); end

%% SMAP with the file
here = pwd; cleanup = onCleanup(@() cd(here));
cd(smapRoot);
SMAP  % starts SMAP, creates g
g.setPar('ROI_restoreparameters', false);
t0 = tic;
loader = plugin('File', 'Load', 'Loader_sml');
loader.attachPar(g.P);
loader.attachLocData(g.locData);
loader.load([], smlFile);
initGuiAfterLoad(loader);
fprintf('loaded %s: %d localizations, %d sites (%.0f s)\n', smlFile, numel(g.locData.loc.xnm), numel(g.locData.SE.sites), toc(t0));
layers = savedLayers(smlFile);
restoreRenderLayers(g, smlFile);
se = g.locData.SE;
se.processors.SEMainGui.make_siteexplorer;
addEvaluator(se, o.evaluator);
if ~isempty(o.excludeSites)
    ids = arrayfun(@(s) s.ID, se.sites);
    for k = find(ismember(ids, o.excludeSites))'
        se.sites(k).annotation.use = false;  % in memory only
    end
    fprintf('left out on request: site(s) %s\n', mat2str(intersect(ids, o.excludeSites)'));
    if ~strcmp(o.selection, 'use'), warning('excludeSites needs ''selection'', ''use'' to take effect'); end
end

%% the plugin
fit = o.fit;
if ischar(fit) && strcmp(fit, 'auto')
    fit = needsFit(se.sites, o.evaluator, o.selection);
end
m = plugin('ROIManager', 'Analyze', 'LocMoFit_fitSortAverage');
m.attachPar(g.P);
m.attachLocData(g.locData);
m.handle = figure('Visible', 'off');
m.makeGui;
p = m.getAllParameters;
p.settingsFile = o.settings; p.loadSettings = 1;
p.evaluator = o.evaluator;
p.doFit = logical(fit);
p.writeRegister = 0;
p.sortPar = o.sortPar;
p.direction.selection = o.direction;
p.sortselection.selection = o.selection;
p.nBins = o.nBins;
p.binBy.selection = o.binBy;
p.binEdges = o.binEdges;
p.cropMode.selection = o.crop;
p.boxSize = o.boxSize;
p.binGap = o.binGap;
p.saveFile = false;  % the script saves its own outputs (below)
p.layout.selection = o.layout;
p.binSpacing = o.binSpacing;
p.removeHotSpots = logical(o.removeHotSpots);
p.hotFactor = o.hotFactor;
p.name = o.name;
if p.doFit
    fprintf('fitting %d sites, about %.0f min\n', numel(se.sites), numel(se.sites) * 16 / 60);
end
res = m.run(p);
if isempty(res), error('LocMoFit_fitSortAverage stopped; see the messages above.'); end
S = res.summary;

%% outputs
writeAverages(g.locData, res.filenumber, outSml, o.name, (height(S) + 1) * res.spacing);
writetable(S, outCsv);
renderAverages(g.locData, res.filenumber, layers, S, res.half, res.spacing, outPng, o);
fprintf('wrote %s\n      %s\n      %s\n', outSml, outPng, outCsv);
out.summary = S;
out.files = {outSml, outPng, outCsv};
out.fitted = p.doFit;
out.spacing = res.spacing;  % bin b is centred at (b * spacing, spacing) nm in the averages file
out.binHalfWidth = res.half;
end

% -------------------------------------------------------------------------

function layers = savedLayers(smlFile)
% Render layers saved in the file: on/off, channels, grouping, active filters.
layers = struct('on', {}, 'channels', {}, 'grouped', {}, 'filters', {});
s = load(smlFile, 'parameters');
if ~isfield(s, 'parameters'), return, end
R = s.parameters.children.guiRender.children;
names = sort(fieldnames(R));
names = names(startsWith(names, 'Layer'));
for k = 1:numel(names)
    n = str2double(names{k}(6:end));
    L = R.(names{k});
    flt = struct('field', {}, 'lo', {}, 'hi', {}, 'invert', {});
    data = L.children.filterTableGui.table.Data;
    for r = 1:size(data, 1)
        if data{r, 7} && ~ismember(data{r, 1}, {'xnm', 'ynm', 'filenumber'})
            inv = size(data, 2) >= 8 && data{r, 8};
            flt(end + 1) = struct('field', data{r, 1}, 'lo', data{r, 2}, 'hi', data{r, 6}, 'invert', inv); %#ok<AGROW>
        end
    end
    onKey = sprintf('layeron%d', n);
    fmt = R.guiFormat;
    layers(n).on = isfield(fmt, onKey) && fmt.(onKey);
    layers(n).channels = double(L.channels(:))';
    layers(n).grouped = logical(L.groupcheck);
    layers(n).filters = flt;
end
end

function restoreRenderLayers(g, smlFile)
% LocMoFit fits on the render layers, so use the file's saved ones.
s = load(smlFile, 'parameters');
if ~isfield(s, 'parameters'), warning('no saved GUI parameters: default render layers'); return, end
try
    g.children.guiRender.setGuiParameters(s.parameters.children.guiRender, true, false);
    g.locData.filter;  % recompute the layer filters from the restored settings
catch e
    warning('run_LocMoFit:layers', 'render layers not restored (%s): fits use the current layers', e.message);
end
end

function addEvaluator(se, name)
ev = se.processors.eval;
for k = 1:4
    if isfield(ev.children, name), return, end
    ev.addmodule('LocMoFitGUI');  % LocMoFitGUI, then LocMoFitGUI_2, ...
end
error('could not add an evaluator named %s', name);
end

function tf = needsFit(sites, name, selection)
% Fit only if none of the sites has a fit; otherwise use the saved fits.
nFit = 0; nUsed = 0;
for k = 1:numel(sites)
    s = sites(k);
    if strcmp(selection, 'use') && ~s.annotation.use, continue, end
    nUsed = nUsed + 1;
    nFit = nFit + (isfield(s.evaluation, name) && isfield(s.evaluation.(name), 'allParsArg'));
end
tf = nFit == 0;
if ~tf && nFit < nUsed
    fprintf('using the saved fits; %d of %d sites have no fit and are left out (fit them, or pass ''fit'', true)\n', nUsed - nFit, nUsed);
end
end

function writeAverages(locData, filenumber, outSml, name, extentNm)
% The averaged file only, as a SMAP _sml.mat (file number 1).
L = locData.loc;
in = L.filenumber == filenumber;
loc = structfun(@(v) v(in), L, 'UniformOutput', false);
loc.filenumber(:) = 1;
file = locData.files.file(filenumber);
file.number = 1;
file.name = [name '_sml.mat'];
for f = {'tif', 'raw', 'image'}
    if isfield(file, f{1}), file.(f{1}) = []; end
end
pixUm = file.info.cam_pixelsize_um;
w = ceil(extentNm / (pixUm(1) * 1000)) + 2;  % field of view covering all bins
file.info.roi = [0 0 w w];
file.info.Width = w;
file.info.Height = w;
saveloc.loc = loc;
saveloc.file = file;
saveloc.history = {};
fileformat.name = 'sml';
save(outSml, 'saveloc', 'fileformat', '-v7.3');
end

function renderAverages(locData, filenumber, layers, S, ~, spacing, outPng, o)
% Gaussian per localization (sigma = locprecnm), each on render layer filters.
% Same layout as the averages file: bin b centred at (b * spacing, spacing),
% with the plugin's gap between neighbours. Vertically, the image is cropped to
% the range of 99.5 % of the localizations plus a margin.
nb = height(S);
px = 2;
if isempty(layers) || ~any([layers.on])
    ch = unique(double(locData.loc.channel(locData.loc.filenumber == filenumber)))';
    layers = arrayfun(@(c) struct('on', true, 'channels', c, 'grouped', true, ...
        'filters', struct('field', {}, 'lo', {}, 'hi', {}, 'invert', {})), ch);
end
layers = layers([layers.on]);
% localizations of each layer, relative to their bin centre
L = struct('dx', {}, 'dy', {}, 'b', {}, 'sig', {});
for k = 1:numel(layers)
    lay = layers(k);
    if lay.grouped, T = locData.grouploc; else, T = locData.loc; end
    sel = T.filenumber == filenumber & ismember(double(T.channel), lay.channels);
    for f = lay.filters
        v = double(T.(f.field));
        inr = v >= f.lo & v <= f.hi;
        if f.invert, inr = ~inr; end
        sel = sel & inr;
    end
    x = double(T.xnm(sel)); y = double(T.ynm(sel));
    b = round(x / spacing);
    L(k).b = b; L(k).dx = x - b * spacing; L(k).dy = y - spacing;
    L(k).sig = max(round(double(T.locprecnm(sel))), 2);
end
r = arrayfun(@(l) abs(l.dy), L, 'UniformOutput', false);
halfY = ceil((prctile(vertcat(r{:}), 99.5) + 30) / 10) * 10;  % nm, with a margin
step = spacing;
xe = (0.5 * step):px:((nb + 0.5) * step);
ye = (step - halfY):px:(step + halfY);
colours = {[0 1 1], 'hot', [0 1 0], [1 0 1]};
rgb = zeros(numel(ye) - 1, numel(xe) - 1, 3);
layerImgs = cell(1, numel(layers));
counts = zeros(1, numel(layers));
for k = 1:numel(layers)
    in = abs(L(k).dx) < step / 2 & abs(L(k).dy) <= halfY;
    counts(k) = numel(in);
    X = L(k).dx(in) + L(k).b(in) * step; Y = L(k).dy(in) + step; sig = L(k).sig(in);
    H = zeros(numel(ye) - 1, numel(xe) - 1);
    for sg = unique(sig)'
        m = sig == sg;
        H = H + imgaussfilt(histcounts2(Y(m), X(m), ye, xe), sg / px);
    end
    c = colours{1 + mod(k - 1, numel(colours))};
    if ischar(c)  % hot: scaled to the maximum, as the colour map is bright already
        I = min(H / max(H(:)), 1);
        map = hot(256);
        layerRgb = reshape(map(1 + round(I(:) * 255), :), [size(I) 3]);
    else
        v = H(H > 1e-3 * max(H(:)));
        I = min(H / prctile(v, 99.5), 1);
        layerRgb = I .* reshape(c, 1, 1, 3);
    end
    layerImgs{k} = layerRgb;
    rgb = rgb + layerRgb;
end
rgb = min(rgb, 1);
desc = arrayfun(@(k) sprintf('layer %d: ch %s, %d locs', k, mat2str(layers(k).channels), counts(k)), 1:numel(layers), 'UniformOutput', false);
titleText = sprintf('%d sites aligned on the fitted base (model 1), sorted by %s %s. %s', ...
    sum(S.n_sites), o.sortPar, o.direction, strjoin(desc, '; '));
drawPanel(rgb, xe, ye, S, step, titleText, outPng, o.sortPar, strcmp(o.binBy, 'distance'));
[d, n] = fileparts(outPng);
for k = 1:numel(layers)
    drawPanel(layerImgs{k}, xe, ye, S, step, [titleText ' (layer ' num2str(k) ' only)'], fullfile(d, sprintf('%s_layer%d.png', n, k)), ...
        o.sortPar, strcmp(o.binBy, 'distance'));
end
end

function drawPanel(rgb, xe, ye, S, spacing, titleText, file, sortPar, byDistance)
nb = height(S);
f = figure('Visible', 'off', 'Color', 'k', 'Position', [0 0 3200 720]);
ax = axes(f, 'Position', [0.005 0.20 0.99 0.66]);
image(ax, xe([1 end]), ye([1 end]), rgb); axis(ax, 'image'); axis(ax, 'off');  % y down, as in SMAP
hold(ax, 'on');
for b = 1:nb
    if byDistance  % the bin's range
        label = sprintf('bin %d (n=%d)\n%s %g..%g nm', b, S.n_sites(b), sortPar, S.range_from(b), S.range_to(b));
    else           % the bin's median
        label = sprintf('bin %d (n=%d)\n%s = %.0f nm', b, S.n_sites(b), sortPar, S.median_sort_value(b));
    end
    text(ax, b * spacing, ye(end) + 8, label, ...
        'Color', 'w', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 18, 'Interpreter', 'none');
end
x0 = (nb + 0.5) * spacing - 140; y0 = ye(1) + 15;
plot(ax, [x0 x0 + 100], [y0 y0], 'w', 'LineWidth', 6);
text(ax, x0 + 50, y0 + 8, '100 nm', 'Color', 'w', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 16);
annotation(f, 'textbox', [0 0.87 1 0.13], 'String', titleText, 'Color', 'w', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 18, 'Interpreter', 'none');
exportgraphics(f, file, 'Resolution', 150, 'BackgroundColor', 'k');
close(f);
end
