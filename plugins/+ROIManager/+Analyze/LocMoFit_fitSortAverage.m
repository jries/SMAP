classdef LocMoFit_fitSortAverage < interfaces.DialogProcessor & interfaces.SEProcessor
% Fit all sites with LocMoFit, sort them by a fitted parameter, align each site
% on its fit and average the sorted sites in a fixed number of bins.
%
% Start from sites that are marked and rotated in the ROI manager, with a
% LocMoFitGUI added in the ROI manager's Evaluate tab. Then, in one run:
%   1. Load the LocMoFit settings file into that LocMoFitGUI (optional).
%   2. Fit all sites ('redraw all' in Evaluate, as clicking through them).
%   3. Sort the sites: usable fits first (use and a successful fit), by the
%      sort parameter (default m2.lPar.y, the base-cap distance), descending by
%      default: flat sites first (bin 1, left), most invaginated last (right).
%   4. Optionally run all sites again with 'Review only' and 'Register',
%      which writes LocMoFit's xnmaligned_/ynmaligned_/rank_ fields in the
%      sorted order. The averages below do not need these fields.
%   5. Align every usable site on its fit and average them in bins, from the
%      first to the last usable site, with the bins as equal as possible.
%
% Alignment: localizations are taken in the site frame exactly as LocMoFit
% fits them (getLocs: subtract site.pos, rotate by the ROI angle), then
% moved by model 1's fitted position, rotation and scale with LocMoFit's own
% locsHandler, as LocMoFit's 'Register' does (run_LocMoFitGUI.m:395-414).
% Model 1's fitted centre (the base) goes to the origin.
%
% Crop ('site ROI', default): each site contributes the localizations of its
% site ROI box around site.pos, exactly those AverageSites_window takes, and
% they are then aligned on the fit, without a second crop. Each site's box lands
% at a slightly different place in the aligned frame, so the edge of the
% average fades out softly instead of showing an outline.
% Crop ('aligned square'): every site is cut to the same square in the aligned
% frame, of side 'box size' (default 1.5 x the site ROI size). This gives every
% site the same support, but a sharp square edge in the averages.
% Bins ('bin by'): 'sites' (default) makes 'bins' bins with equal numbers of
% sites, from the first to the last usable site. 'distance' makes one bin per
% range of the sort value between consecutive 'bin edges' (nm), e.g.
% Inf,15,0,-10,-20,-35,-50,-70,-95,-125,-Inf; the ranges are then the same for
% every data set, bins appear left to right in the order of the edges, and sites
% outside all ranges are left out. The csv lists each bin's range.
%
% Layout ('fixed', default): the fitted base (model 1) of every site is placed at
% (b x spacing, spacing) for bin b, with a fixed spacing (default 550 nm). The
% base is then at the same coordinates in every data set, so the averages of
% different data sets (e.g. different proteins against the same reference) can
% be overlaid directly. Use the same spacing for all of them.
% Layout ('compact'): the aligned frame is shifted so that the ROI centres lie on
% average at the bin centre, and the bins are placed with 'gap' nm (default
% 250) between neighbours. Spacing and centring then depend on the data set.
%
% Hot spots (on by default): a place where one site has far more localizations
% within 12 nm than sites normally have at their densest place (a dye stuck on,
% a bead, dirt) shows up as a bright dot in the average. Per channel, the
% threshold is 'factor' (default 3) x the 90th percentile over all sites of
% their densest count, and at least 50; that site's localizations within 20 nm
% of the spot are left out of the average, the rest of the site stays. The
% removed sites are listed in the command window.
%
% Output: one new file with the bins side by side (bin 1 left), class = site
% ID, and a summary in the command window and the results tab. With 'save to
% source folder' (default on), the averages are also saved, at the end of the
% run, as <source>_aligned_avg_<n>bins_sml.mat with the bin table as .csv and
% the sort value against the sorted sites as <name>_<sortPar>_vs_site.png (with
% its points as .csv), in the folder of the loaded _sml.mat (lastSMLFile);
% nothing is overwritten.
%
% Default LocMoFit settings: LocMoFit_settings/Myo5_Sla2_settings_ver2_LocMoFit.mat
% next to this file (a copy of the settings used for Myo5/Sla2; the models are
% in LocMoFit/models/obsolete). They are printed in the command window when
% loaded. Optimizer particleswarm (Display iter, UseVectorized 1).
%   model 1: thickRing2DPro_continuous (base), layer 1, pixel 5 nm, sigma 15 nm
%   model 2: hemiellipse2DPro_continuous (cap, relative to model 1), layer 1,
%            pixel 5 nm, sigma 15 nm
%   lb/ub are relative to the start value, min/max absolute:
%   model type   name         fixed    value     lb     ub    min    max
%   1     lPar   x            no      -21.37    -30     30   -150    150
%   1     lPar   y            no       79.45      0    125   -150    150
%   1     lPar   zrot         no       -4.78    -25     25   -Inf    Inf
%   1     lPar   variation    yes          0      0     10      0     20
%   1     lPar   x/yscale     yes          1      1      1      1      1
%   1     lPar   weight       yes          1  1e-05      1  1e-05      1
%   1     mPar   xcenter      yes          0    -60     60    -60     60
%   1     mPar   ycenter      yes          0     30    230     30    230
%   1     mPar   innerRadius  yes         20    -10     10     20     80
%   1     mPar   outerRadius  yes         90    -10     10     80    100
%   1     mPar   thickness    yes         60     -5      5     40     70
%   2     lPar   x            yes          0   -150    150   -150    150
%   2     lPar   y            no      -17.44   -260     20   -300    200
%   2     lPar   zrot         yes          0   -Inf    Inf   -Inf    Inf
%   2     lPar   variation    yes          0      0     10      0     20
%   2     lPar   x/yscale     yes          1      1      1      1      1
%   2     lPar   weight       yes          1  1e-05      1  1e-05      1
%   2     mPar   xcenter      yes          0    -60     60   -150    150
%   2     mPar   ycenter      yes          0   -230     30   -150    150
%   2     mPar   a            yes         40    -10     10     10     40
%   2     mPar   b            no         250      0     40      0    270
%   91    offset weight       no      0.6115     -1      1  0.001  0.999
%   Converter: start value of m1.lPar.y = median(locs.ynm(locs.layer==1)) per site.
    methods
        function obj = LocMoFit_fitSortAverage(varargin)
            obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters = {'se_siteroi'};
            obj.showresults = true;
        end

        function out = run(obj, p)
            out = [];
            global SMAP_stopnow
            se = obj.SE;
            if isempty(se.sites)
                errordlg('No sites in the ROI manager.', 'LocMoFit_fitSortAverage'); return
            end
            if ~isfield(se.processors, 'preview') || ~isvalid(se.processors.preview)
                se.processors.SEMainGui.make_siteexplorer;  % the site viewer, needed by 'redraw all'
            end
            evalName = strtrim(p.evaluator);
            gui = findLocMoFitGUI(se, evalName);
            if isempty(gui), return, end

            % 1. settings
            if p.loadSettings
                file = strtrim(p.settingsFile);
                if ~isfile(file)
                    errordlg(['Settings file not found: ' file], 'LocMoFit_fitSortAverage'); return
                end
                obj.status('loading LocMoFit settings'); drawnow
                gui.load_callback([], [], file);
                printFitSettings(file, gui.fitter);
            end
            if isempty(gui.fitter)
                errordlg(['No model loaded in ' evalName '. Load the settings first.'], 'LocMoFit_fitSortAverage'); return
            end

            % 2. fit all sites
            if p.doFit
                warnOtherEvaluators(se, evalName);
                setFitMode(gui, false, false);
                nsites = numel(se.sites);
                selected = se.processors.preview.guihandles.sitelist.Value;
                if numel(selected) > 1, nsites = numel(selected); end  % redraw all then fits only these
                obj.status(sprintf('fitting %d sites (Evaluate: redraw all)', nsites)); drawnow
                t0 = tic;
                se.processors.eval.redrawall;
                fprintf('LocMoFit_fitSortAverage: fitted %d sites in %.1f min\n', nsites, toc(t0)/60);
                if SMAP_stopnow
                    obj.status('stopped by the user after fitting'); return
                end
            end

            % 3. sort
            sites = se.sites;
            [usable, why] = fitStatus(sites, evalName, p.sortselection.selection);
            keyIdx = parameterIndex(sites(find(usable, 1)), evalName, p.sortPar);
            if isempty(keyIdx)
                errordlg(['Sort parameter not found: ' p.sortPar], 'LocMoFit_fitSortAverage'); return
            end
            key = nan(numel(sites), 1);
            for k = find(usable)'
                key(k) = sites(k).evaluation.(evalName).allParsArg.value(keyIdx);
            end
            sortkey = key;
            if strcmp(p.direction.selection, 'descend'), sortkey = -sortkey; end
            sortkey(~usable) = inf;
            [~, order] = sortrows([double(~usable(:)) sortkey]);
            se.sites = sites(order);
            se.setIndList;
            se.processors.preview.updateSitelist;
            sites = se.sites;
            usable = usable(order);
            key = key(order);
            n = sum(usable);  % the usable sites are now 1..n
            reportExcluded(why(order), numel(sites), n);
            % which bin each usable site goes to (0: none)
            [binOfSite, binRange] = assignBins(key(1:n), p);
            if isempty(binRange), return, end
            nB = size(binRange, 1);

            % 4. optional: SMAP's own Register fields, in the sorted order
            if p.writeRegister
                setFitMode(gui, true, true);
                obj.status('Review only + Register for all sites'); drawnow
                se.processors.eval.redrawall;
                setFitMode(gui, false, false);
                if SMAP_stopnow
                    obj.status('stopped by the user during Register'); return
                end
            end

            % 5. align and average in bins
            boxSize = p.boxSize;
            if isempty(boxSize) || ~isfinite(boxSize) || boxSize <= 0
                boxSize = 1.5 * p.se_siteroi;  % room for long invaginations above the base
            end
            fitter = gui.fitter;
            tf = siteTransforms(sites(1:n), evalName, fitter);
            fixedLayout = ~isfield(p, 'layout') || strcmp(p.layout.selection, 'fixed');
            if fixedLayout
                centre = [0 0];  % the fitted base (model 1) stays at the bin centre
            else
                centre = median(tf.roiCentre, 1);  % where the ROI centres lie, aligned frame
            end
            roiCrop = strcmp(p.cropMode.selection, 'site ROI');
            if roiCrop
                % each site's localizations as AverageSites_window takes them: the
                % site ROI box around site.pos (getLocs, 'size', [se_siteroi se_siteroi])
                region = @(k) [sites(k).pos(1) sites(k).pos(2) p.se_siteroi p.se_siteroi];
            else
                % a circle around site.pos that contains the aligned square of every site
                half = boxSize / 2;
                readRadius = tf.maxScale * (hypot(centre(1), centre(2)) + half * sqrt(2)) + tf.shift + 1;
                region = @(k) [sites(k).pos(1) sites(k).pos(2) readRadius(k)];
            end

            newfile = obj.locData.files.filenumberEnd + 1;
            fileName = sprintf('%s_%dbins_F%d', strtrim(p.name), nB, newfile);
            obj.locData.addfile(fileName);
            initGuiAfterLoad(obj);
            se.processors.preview.updateFilelist;

            loc = obj.locData.loc;
            if ~isfield(loc, 'class')
                obj.locData.setloc('class', zeros(size(loc.xnm), 'single'));
                loc = obj.locData.loc;
            end
            fields = fieldnames(loc);
            % bin, first site, last site, sites, median key, locs, range from, range to
            summary = zeros(nB, 8);
            bins = struct('ind', cell(nB, 1), 'x', [], 'y', [], 'cls', []);
            ticc = tic;
            for b = 1:nB
                ks = find(binOfSite == b)';
                ind = cell(numel(ks), 1); xa = ind; ya = ind; cls = ind;
                for j = 1:numel(ks)
                    k = ks(j);
                    [ind{j}, r] = alignedSiteLocs(obj.locData, sites(k), tf, k, region(k), fitter);
                    if roiCrop
                        keep = true(size(ind{j}));
                    else
                        keep = abs(r(:, 1) - centre(1)) <= half & abs(r(:, 2) - centre(2)) <= half;
                    end
                    ind{j} = ind{j}(keep);
                    xa{j} = r(keep, 1) - centre(1);
                    ya{j} = r(keep, 2) - centre(2);
                    cls{j} = repmat(sites(k).ID, nnz(keep), 1);
                    if toc(ticc) > 1
                        obj.status(sprintf('averaging: bin %d of %d, site %d of %d', b, nB, k, n)); drawnow
                        ticc = tic;
                    end
                end
                bins(b).ind = vertcat(ind{:});
                bins(b).x = vertcat(xa{:});
                bins(b).y = vertcat(ya{:});
                bins(b).cls = vertcat(cls{:});
                if isempty(ks)
                    summary(b, :) = [b, NaN, NaN, 0, NaN, 0, binRange(b, :)];
                else
                    summary(b, :) = [b, ks(1), ks(end), numel(ks), median(key(ks)), numel(bins(b).ind), binRange(b, :)];
                end
            end
            % hot spots: one site emitting again and again at one place (a dye stuck
            % on, a bead, dirt) shows up as a bright dot in the average
            if p.removeHotSpots
                binOf = repelem((1:nB)', arrayfun(@(bn) numel(bn.ind), bins));
                allInd = vertcat(bins.ind);
                keepAll = hotSpotKeep(vertcat(bins.x), vertcat(bins.y), vertcat(bins.cls), ...
                    double(loc.channel(allInd)), p.hotFactor);
                for b = 1:nB
                    k = keepAll(binOf == b);
                    bins(b).ind = bins(b).ind(k); bins(b).x = bins(b).x(k);
                    bins(b).y = bins(b).y(k); bins(b).cls = bins(b).cls(k);
                end
            end
            % bins side by side
            if roiCrop
                % where the averages visibly end: 95 % of the localizations
                half = prctile(abs(vertcat(bins.x)), 95);
            end
            if fixedLayout
                % a fixed spacing, the same for every data set, so that averages of
                % different data sets can be overlaid: base of bin b at (b*S, S)
                spacing = p.binSpacing;
                if isempty(spacing) || ~isfinite(spacing) || spacing <= 0, spacing = 550; end
            else
                % binGap nm of empty space between neighbours (depends on the data)
                gap = p.binGap;
                if isempty(gap) || ~isfinite(gap) || gap < 0, gap = 250; end
                spacing = 2 * half + gap;
            end
            nOut = sum(arrayfun(@(bn) nnz(abs(bn.x) >= spacing / 2), bins));
            if nOut > 0
                fprintf('LocMoFit_fitSortAverage: %d of %d localizations (%.2f %%) are further than %.0f nm from their bin centre and left out.\n', ...
                    nOut, numel(vertcat(bins.x)), 100 * nOut / numel(vertcat(bins.x)), spacing / 2);
            end
            for b = 1:nB
                % within half a spacing of the bin centre, so bins never mix
                in = abs(bins(b).x) < spacing / 2;
                summary(b, 6) = nnz(in);
                locnew = struct();
                for f = 1:numel(fields)
                    locnew.(fields{f}) = loc.(fields{f})(bins(b).ind(in), :);
                end
                locnew.xnm = bins(b).x(in) + b * spacing;
                locnew.ynm = bins(b).y(in) + spacing;
                locnew.class = bins(b).cls(in);
                locnew.filenumber = repmat(newfile, nnz(in), 1);
                for f = 1:numel(fields)
                    obj.locData.addloc(fields{f}, locnew.(fields{f}));
                end
            end
            obj.locData.regroup;
            obj.locData.filter;

            printSummary(summary, p.sortPar, fileName, 2 * half, centre);
            plotSummary(obj.initaxis('bins'), key(1:n), binOfSite, binRange, p.sortPar);
            obj.status(sprintf('done: %d bins of %d sites in file %d (%s)', nB, nnz(binOfSite), newfile, fileName));
            out.summary = array2table(summary, 'VariableNames', ...
                {'bin', 'first_site', 'last_site', 'n_sites', 'median_sort_value', 'n_locs', 'range_from', 'range_to'});
            out.filenumber = newfile;
            out.spacing = spacing;  % bin b is centred at (b * spacing, spacing)
            out.centre = centre;    % the fitted base lies at (b * spacing, spacing) - centre
            out.half = half;        % half-width (nm) of each average (95 % of its localizations)

            % 6. save the averages next to the source file
            if ~isfield(p, 'saveFile') || p.saveFile
                byDistance = isfield(p, 'binBy') && strcmp(p.binBy.selection, 'distance');
                curve = struct('key', key(1:n), 'bin', binOfSite, 'range', binRange, 'sortPar', p.sortPar, ...
                    'id', arrayfun(@(st) st.ID, sites(1:n)));
                out.savedFiles = saveAverages(obj, newfile, out.summary, nB, spacing, byDistance, curve);
            end
        end

        function pard = guidef(obj)
            pard = guidef(obj);
        end
    end
end

% -------------------------------------------------------------------------
% Steps
% -------------------------------------------------------------------------

function gui = findLocMoFitGUI(se, evalName)
gui = [];
children = se.processors.eval.children;
if isfield(children, evalName)
    gui = children.(evalName);
    return
end
names = fieldnames(children);
names = names(startsWith(names, 'LocMoFit'));
if isempty(names)
    msg = 'Add a LocMoFitGUI in the ROI manager''s Evaluate tab first.';
else
    msg = sprintf('No evaluator named %s. Loaded: %s', evalName, strjoin(names, ', '));
end
errordlg(msg, 'LocMoFit_fitSortAverage');
end

function warnOtherEvaluators(se, evalName)
% 'redraw all' runs every switched-on evaluator, so other ones cost time too.
data = se.processors.eval.guihandles.modules.Data;
if size(data, 2) < 2, return, end
on = data([data{:, 1}] == true, 2);
others = setdiff(on, {evalName});
if ~isempty(others)
    warning('LocMoFit_fitSortAverage:otherEvaluators', ...
        'These evaluators also run for every site: %s', strjoin(others, ', '));
end
selected = se.processors.preview.guihandles.sitelist.Value;
if numel(selected) > 1
    warning('LocMoFit_fitSortAverage:selection', ...
        '%d sites are selected in the site list: redraw all fits only those.', numel(selected));
end
end

function setFitMode(gui, reviewOnly, register)
% The 'Review only' (noFit) and 'Register' (useAlignment) checkboxes of LocMoFitGUI.
gui.guihandles.noFit.Value = reviewOnly;
gui.guihandles.useAlignment.Value = register;
end

function [usable, why] = fitStatus(sites, evalName, selection)
% A site is usable if it has a successful fit and, for selection 'use', use checked.
usable = false(numel(sites), 1);
why = repmat({''}, numel(sites), 1);
for k = 1:numel(sites)
    s = sites(k);
    if strcmp(selection, 'use') && ~s.annotation.use
        why{k} = 'use unchecked';
    elseif ~isfield(s.evaluation, evalName) || ~isfield(s.evaluation.(evalName), 'allParsArg')
        why{k} = 'not fitted';
    else
        ev = s.evaluation.(evalName);
        a = ev.allParsArg;
        if ~isstruct(a) || ~isfield(a, 'value') || ~isfield(a, 'fix') || isempty(a.value) || numel(a.value) ~= numel(a.fix)
            why{k} = 'incomplete fit parameters';  % e.g. a fit that stopped before saving its values
            continue
        end
        free = ~a.fix;
        if isfield(ev, 'fitInfo') && isfield(ev.fitInfo, 'guiInfo') && strcmp(ev.fitInfo.guiInfo, 'Fit or plot failed.')
            why{k} = 'fit failed';
        elseif ~all(isfinite(ev.allParsArg.value(free)))
            why{k} = 'non-finite fit parameters';
        else
            usable(k) = true;
        end
    end
end
end

function idx = parameterIndex(site, evalName, spec)
% Row of allParsArg for 'm2.lPar.y', 'm2.y', 'value(14)' or '14'.
idx = [];
if isempty(site), return, end
spec = strtrim(spec);
num = str2double(regexprep(spec, '^value\((\d+)\)$', '$1'));
a = site.evaluation.(evalName).allParsArg;
if isfinite(num)
    if num >= 1 && num <= numel(a.value), idx = num; end
    return
end
parts = strsplit(spec, '.');  % m2.lPar.y or m2.y
if numel(parts) < 2 || numel(parts) > 3 || isempty(regexp(parts{1}, '^m\d+$', 'once')), return, end
hit = a.model == str2double(parts{1}(2:end)) & strcmp(a.name, parts{end});
if numel(parts) == 3, hit = hit & strcmp(a.type, parts{2}); end
idx = find(hit, 1);
end

function tf = siteTransforms(sites, evalName, fitter)
% Model 1's location parameters per site, and where each ROI centre lands
% in the aligned frame (the registration of the site-frame origin).
n = numel(sites);
tf.lPar = cell(n, 1);
tf.roiCentre = zeros(n, 2);
tf.shift = zeros(n, 1);
tf.maxScale = zeros(n, 1);
for k = 1:n
    a = sites(k).evaluation.(evalName).allParsArg;
    lp = struct();
    rows = find(a.model == 1 & strcmp(a.type, 'lPar'));
    for r = rows(:)'
        lp.(a.name{r}) = a.value(r);
    end
    tf.lPar{k} = lp;
    o = fitter.locsHandler(struct('xnm', 0, 'ynm', 0), lp, 1);
    tf.roiCentre(k, :) = [o.xnm, o.ynm];
    tf.shift(k) = hypot(lp.x, lp.y);
    tf.maxScale(k) = max([getfield_default(lp, 'xscale', 1), getfield_default(lp, 'yscale', 1)]);
end
end

function v = getfield_default(s, name, default)
if isfield(s, name), v = s.(name); else, v = default; end
end

function [ind, r] = alignedSiteLocs(locData, site, tf, k, region, fitter)
% Localizations of the site's file in region ([x y w h]: box of width w and
% height h, [x y r]: circle of radius r), ungrouped and unfiltered, in model
% 1's aligned frame. ind indexes locData.loc.
[locs, inside] = locData.getloc({'xnm', 'ynm', 'filenumber'}, ...
    'position', region, 'grouping', 'ungrouped');
ind = find(inside);
own = locs.filenumber == site.info.filenumber;
ind = ind(own);
% site frame, as SEEvaluationProcessor.getLocs (xnmrot, ynmrot)
[q.xnm, q.ynm] = rotcoord(double(locs.xnm(own)) - site.pos(1), double(locs.ynm(own)) - site.pos(2), ...
    site.annotation.rotationpos.angle / 180 * pi);
% model 1's frame, as LocMoFit's Register (locsRegister -> locsHandler)
a = fitter.locsHandler(q, tf.lPar{k}, 1);
r = [a.xnm(:), a.ynm(:)];
end

function keep = hotSpotKeep(x, y, site, ch, factor)
% Remove hot spots: places where one site has far more localizations within
% 12 nm than sites normally have at their densest place. Per channel, the
% threshold is factor x the 90th percentile over sites of that densest count,
% and at least 50. Localizations of the site within 20 nm of a hot spot go.
rCount = 12; rRemove = 20; px = 3;  % nm
if isempty(factor) || ~isfinite(factor) || factor <= 0, factor = 3; end
keep = true(size(x));
count = zeros(size(x));
[g, gSite, gCh] = findgroups(site, ch);
peak = zeros(max(g), 1);
for k = 1:max(g)
    m = find(g == k);
    count(m) = neighbourCount(x(m), y(m), rCount, px);
    peak(k) = max(count(m));
end
removed = zeros(0, 4);  % site, channel, densest count, removed locs
for c = unique(gCh)'
    thr = max(50, factor * prctile(peak(gCh == c), 90));
    for k = find(gCh == c & peak > thr)'
        m = find(g == k);
        hot = count(m) > thr;
        [~, near] = neighbourCount(x(m), y(m), rRemove, px, hot);
        keep(m(near)) = false;
        removed(end + 1, :) = [gSite(k), c, peak(k), nnz(near)]; %#ok<AGROW>
    end
    fprintf('hot spots, channel %g: threshold %.0f localizations within %d nm (densest spot per site: median %.0f, 90th pct %.0f)\n', ...
        c, thr, rCount, median(peak(gCh == c)), prctile(peak(gCh == c), 90));
end
if isempty(removed)
    fprintf('hot spots: none found\n');
else
    removed = sortrows(removed, -4);
    fprintf('hot spots removed (%d localizations, %.2f %%):\n  site  channel  densest  removed\n', sum(removed(:, 4)), 100 * sum(removed(:, 4)) / numel(x));
    fprintf('  %4d  %7g  %7d  %7d\n', removed');
end
end

function [count, near] = neighbourCount(x, y, r, px, hot)
% count: localizations within about r of each localization (on a px grid).
% near: with hot given, the localizations within about r of any hot one.
ix = floor((x - min(x)) / px) + 1;
iy = floor((y - min(y)) / px) + 1;
sz = [max(ix), max(iy)];
n = ceil(r / px);
[dx, dy] = meshgrid(-n:n);
disk = double(dx.^2 + dy.^2 <= n^2);
lin = sub2ind(sz, ix, iy);
if nargin < 5
    H = accumarray(lin, 1, [prod(sz), 1]);
    C = conv2(reshape(H, sz), disk, 'same');
    count = C(lin);
    near = [];
else
    H = accumarray(lin(hot), 1, [prod(sz), 1]);
    D = conv2(reshape(H, sz), disk, 'same') > 0;
    near = D(lin);
    count = [];
end
end

% -------------------------------------------------------------------------
% Reporting
% -------------------------------------------------------------------------

function reportExcluded(why, nAll, nUsable)
fprintf('LocMoFit_fitSortAverage: %d of %d sites used', nUsable, nAll);
reasons = why(~cellfun(@isempty, why));
if isempty(reasons)
    fprintf('.\n'); return
end
[u, ~, j] = unique(reasons);
counts = accumarray(j, 1);
parts = arrayfun(@(i) sprintf('%d %s', counts(i), u{i}), 1:numel(u), 'UniformOutput', false);
fprintf('; left out (moved to the end of the list): %s.\n', strjoin(parts, ', '));
end

function saved = saveAverages(obj, filenumber, summary, nB, spacing, byDistance, curve)
% The averages file entry as a SMAP _sml.mat, the bin table as .csv, and the
% sort value against the sorted sites (the 'bins' plot) as .png with its points
% as .csv, in the folder of the last loaded _sml.mat; a name that exists gets
% _2, _3, ...
saved = {};
src = obj.getPar('lastSMLFile');
if isempty(src) || ~isfile(src)
    warning('LocMoFit_fitSortAverage:noSource', 'no loaded _sml.mat known: the averages were not saved');
    return
end
[folder, stem] = fileparts(src);
base = sprintf('%s_aligned_avg_%dbins', regexprep(stem, '_sml$', ''), nB);
if byDistance, base = [base '_bydistance']; end
name = base;
k = 1;
tag = ['_' regexprep(curve.sortPar, '[^A-Za-z0-9]', '_') '_vs_site'];   % e.g. _m2_lPar_y_vs_site
while isfile(fullfile(folder, [name '_sml.mat'])) || isfile(fullfile(folder, [name '.csv'])) ...
        || isfile(fullfile(folder, [name tag '.png']))
    k = k + 1;
    name = sprintf('%s_%d', base, k);
end
outSml = fullfile(folder, [name '_sml.mat']);
outCsv = fullfile(folder, [name '.csv']);
obj.status('saving the averages'); drawnow
L = obj.locData.loc;
in = L.filenumber == filenumber;
saveloc.loc = structfun(@(v) v(in), L, 'UniformOutput', false);
saveloc.loc.filenumber(:) = 1;
file = obj.locData.files.file(filenumber);
file.number = 1;
file.name = [name '_sml.mat'];
for f = {'tif', 'raw', 'image'}
    if isfield(file, f{1}), file.(f{1}) = []; end
end
w = ceil((nB + 1) * spacing / (file.info.cam_pixelsize_um(1) * 1000)) + 2;  % field of view over all bins
file.info.roi = [0 0 w w];
file.info.Width = w;
file.info.Height = w;
saveloc.file = file;
saveloc.history = {};
fileformat.name = 'sml';
save(outSml, 'saveloc', 'fileformat', '-v7.3');
writetable(summary, outCsv);
outPng = fullfile(folder, [name tag '.png']);
outCurve = fullfile(folder, [name tag '.csv']);
f = figure('Visible', 'off', 'Position', [100 100 900 500]);
ax = axes(f);
plotSummary(ax, curve.key, curve.bin, curve.range, curve.sortPar);
title(ax, sprintf('%s: %d sites in %d bins', strrep(name, '_', '\_'), numel(curve.key), nB));
exportgraphics(f, outPng, 'Resolution', 200);
close(f);
writetable(table((1:numel(curve.key))', curve.id(:), curve.key(:), curve.bin(:), ...
    'VariableNames', {'site_rank', 'site_id', regexprep(curve.sortPar, '[^A-Za-z0-9]', '_'), 'bin'}), outCurve);
saved = {outSml, outCsv, outPng, outCurve};
fprintf('saved %s\n      %s\n      %s\n      %s\n', outSml, outCsv, outPng, outCurve);
end

function printSummary(summary, sortPar, fileName, extent, centre)
fprintf('Averages in file %s (each bin within %.0f nm; base at (%.0f, %.0f) nm from the bin centre):\n', ...
    fileName, extent, -centre(1), -centre(2));
fprintf('  bin   sites      n   median %-10s  locs      %s range\n', sortPar, sortPar);
for b = 1:size(summary, 1)
    fprintf('  %3d   %4d-%-4d %4d   %10.1f   %9d   %g .. %g\n', summary(b, :));
end
end

function plotSummary(ax, key, binOfSite, binRange, sortPar)
plot(ax, 1:numel(key), key, '.');
hold(ax, 'on');
for e = find(diff(binOfSite(:)) ~= 0)'  % where the bin changes along the sorted sites
    xline(ax, e + 0.5, ':');
end
hold(ax, 'off');
xlabel(ax, 'site (sorted)');
ylabel(ax, sortPar);
title(ax, sprintf('%d bins', size(binRange, 1)));
end

function [binOfSite, binRange] = assignBins(key, p)
% Bin of each usable site (sorted order), and each bin's range of the sort value.
% 'sites': nBins bins with (nearly) equal numbers of sites, from the first to the
% last usable site. 'distance': bins between consecutive binEdges (sort value),
% lower edge included, in the order the edges are given (left to right); sites
% outside all bins are left out.
n = numel(key);
binOfSite = zeros(n, 1);
binRange = zeros(0, 2);
byDistance = isfield(p, 'binBy') && strcmp(p.binBy.selection, 'distance');
if ~byDistance
    if n < p.nBins
        errordlg(sprintf('Only %d usable sites for %d bins.', n, p.nBins), 'LocMoFit_fitSortAverage'); return
    end
    e = round(linspace(0, n, p.nBins + 1));
    binRange = zeros(p.nBins, 2);
    for b = 1:p.nBins
        binOfSite(e(b) + 1:e(b + 1)) = b;
        binRange(b, :) = [key(e(b) + 1), key(e(b + 1))];  % first and last site's value
    end
    return
end
edges = p.binEdges;
if ischar(edges) || isstring(edges), edges = str2num(['[' char(edges) ']']); end %#ok<ST2NM>
edges = double(edges(:))';
d = diff(edges);
if numel(edges) < 2 || any(isnan(edges)) || ~(all(d > 0) || all(d < 0))
    errordlg('bin edges: give at least two values, all increasing or all decreasing, e.g. Inf,15,0,-10,-125,-Inf', 'LocMoFit_fitSortAverage');
    return
end
nB = numel(edges) - 1;
binRange = [edges(1:end-1)', edges(2:end)'];
for b = 1:nB
    lo = min(binRange(b, :)); hi = max(binRange(b, :));
    in = key >= lo & key < hi;
    if hi == Inf, in = key >= lo; end
    binOfSite(in & binOfSite == 0) = b;
end
nOut = nnz(binOfSite == 0);
counts = accumarray(binOfSite(binOfSite > 0), 1, [nB, 1])';
fprintf('bins by %s: edges %s; sites per bin %s', p.sortPar, mat2str(edges), mat2str(counts));
if nOut > 0, fprintf('; %d site(s) outside all bins, left out', nOut); end
fprintf('\n');
if any(counts == 0)
    warning('LocMoFit_fitSortAverage:emptyBin', 'bin(s) %s have no sites', mat2str(find(counts == 0)));
end
end

function browse_callback(~, ~, obj)
f = obj.getSingleGuiParameter('settingsFile');
[file, path] = uigetfile({'*_LocMoFit.mat;*_SMLMModelFit.mat', 'LocMoFit settings'}, 'LocMoFit settings', f);
if isequal(file, 0), return, end
obj.setGuiParameters(struct('settingsFile', fullfile(path, file)));
end

% -------------------------------------------------------------------------
% GUI
% -------------------------------------------------------------------------

function f = defaultSettingsFile
% the settings that come with the plugin (LocMoFit_settings, next to this file)
f = fullfile(fileparts(mfilename('fullpath')), 'LocMoFit_settings', 'Myo5_Sla2_settings_ver2_LocMoFit.mat');
end

function printFitSettings(file, fitter)
% What will be fitted with: models and parameter table of the loaded settings.
fprintf('LocMoFit settings: %s\n  optimizer %s\n', file, fitter.solver.SolverName);
for m = 1:fitter.numOfModel
    fprintf('  model %d: %s, layer %d\n', m, class(fitter.model{m}.modelObj), fitter.model{m}.layer);
end
a = fitter.allParsArg;
fprintf('  %-6s %-7s %-12s %-6s %9s %8s %8s %8s %8s\n', 'model', 'type', 'name', 'fixed', 'value', 'lb', 'ub', 'min', 'max');
yn = {'no', 'yes'};
for r = 1:numel(a.value)
    fprintf('  %-6d %-7s %-12s %-6s %9.4g %8.4g %8.4g %8.4g %8.4g\n', a.model(r), a.type{r}, a.name{r}, ...
        yn{1 + a.fix(r)}, a.value(r), a.lb(r), a.ub(r), a.min(r), a.max(r));
end
end

function pard = guidef(obj)
pard.loadSettings.object = struct('String', 'load settings', 'Style', 'checkbox', 'Value', 1);
pard.loadSettings.position = [1, 1];
pard.loadSettings.Width = 1;
pard.settingsFile.object = struct('String', defaultSettingsFile, 'Style', 'edit');
pard.settingsFile.position = [1, 2];
pard.settingsFile.Width = 2.4;
pard.settingsFile.TooltipString = 'LocMoFit settings (*_LocMoFit.mat), saved with the save button of LocMoFitGUI. Default: the Myo5/Sla2 settings in LocMoFit_settings next to this plugin.';
pard.browse.object = struct('String', '...', 'Style', 'pushbutton', 'Callback', {{@browse_callback, obj}});
pard.browse.position = [1, 4.4];
pard.browse.Width = 0.6;

pard.tEvaluator.object = struct('String', 'evaluator', 'Style', 'text');
pard.tEvaluator.position = [2, 1];
pard.evaluator.object = struct('String', 'LocMoFitGUI', 'Style', 'edit');
pard.evaluator.position = [2, 2];
pard.evaluator.TooltipString = 'Name of the LocMoFitGUI in the Evaluate tab (e.g. LocMoFitGUI or LocMoFitGUI_2).';
pard.doFit.object = struct('String', 'fit all sites', 'Style', 'checkbox', 'Value', 1);
pard.doFit.position = [2, 3];
pard.doFit.Width = 2;
pard.doFit.TooltipString = 'Fit every site (Evaluate: redraw all). Untick to use the fits already in the sites.';

pard.tSort.object = struct('String', 'sort by', 'Style', 'text');
pard.tSort.position = [3, 1];
pard.sortPar.object = struct('String', 'm2.lPar.y', 'Style', 'edit');
pard.sortPar.position = [3, 2];
pard.sortPar.TooltipString = 'Fit parameter: m2.lPar.y (base-cap distance), m2.y, value(14) or 14.';
pard.direction.object = struct('String', {{'descend', 'ascend'}}, 'Style', 'popupmenu');
pard.direction.position = [3, 3];
pard.direction.TooltipString = ['descend: largest value first. For m2.lPar.y (0 or above for flat sites, negative when ' ...
    'invaginated) this puts flat sites in bin 1 (left) and the most invaginated in the last bin (right).'];
pard.sortselection.object = struct('String', {{'use', 'all'}}, 'Style', 'popupmenu');
pard.sortselection.position = [3, 4];
pard.sortselection.TooltipString = 'use: only sites with use checked; all: every site with a successful fit.';

pard.tBins.object = struct('String', 'bins', 'Style', 'text');
pard.tBins.position = [4, 1];
pard.nBins.object = struct('String', '10', 'Style', 'edit');
pard.nBins.TooltipString = 'bin by sites: number of bins (equal numbers of sites).';
pard.nBins.position = [4, 2];
pard.nBins.Width = 0.5;
pard.cropMode.object = struct('String', {{'site ROI', 'aligned square'}}, 'Style', 'popupmenu');
pard.cropMode.position = [4, 2.6];
pard.cropMode.Width = 1.2;
pard.cropMode.TooltipString = ['site ROI: each site contributes the localizations of its site ROI box, as in AverageSites_window, ' ...
    'then aligned on the fit; the edge fades out softly. aligned square: every site cut to the same square after ' ...
    'alignment (sharp square edge).'];
pard.boxSize.object = struct('String', '', 'Style', 'edit');
pard.boxSize.position = [4, 3.8];
pard.boxSize.Width = 0.5;
pard.boxSize.TooltipString = 'aligned square only: side (nm) of the square. Empty: 1.5 x the site ROI size (se_siteroi).';
pard.removeHotSpots.object = struct('String', 'remove hot spots, factor', 'Style', 'checkbox', 'Value', 1);
pard.removeHotSpots.position = [6, 3];
pard.removeHotSpots.Width = 1.5;
pard.removeHotSpots.TooltipString = ['Remove places where one site emits again and again (dye stuck on, bead, dirt): more localizations within 12 nm ' ...
    'than factor x the 90th percentile over sites of their densest spot (per channel, at least 50). Removed: that site''s localizations within 20 nm.'];
pard.hotFactor.object = struct('String', '3', 'Style', 'edit');
pard.hotFactor.position = [6, 4.5];
pard.hotFactor.Width = 0.4;
pard.saveFile.object = struct('String', 'save to source folder', 'Style', 'checkbox', 'Value', 1);
pard.saveFile.position = [9, 1];
pard.saveFile.Width = 2;
pard.saveFile.TooltipString = ['At the end, save the averages as <source>_aligned_avg_<n>bins_sml.mat, the bin table as .csv and ' ...
    'the sort value against the sorted sites as _vs_site.png and .csv, in the folder of the loaded _sml.mat; an existing name gets _2, _3, ...'];
pard.tBinBy.object = struct('String', 'bin by', 'Style', 'text');
pard.tBinBy.position = [8, 1];
pard.binBy.object = struct('String', {{'sites', 'distance'}}, 'Style', 'popupmenu');
pard.binBy.position = [8, 2];
pard.binBy.Width = 0.8;
pard.binBy.TooltipString = ['sites: ''bins'' bins with equal numbers of sites. distance: bins between the bin edges ' ...
    '(values of the sort parameter, nm), the same ranges for every data set.'];
pard.binEdges.object = struct('String', 'Inf,15,0,-10,-20,-35,-50,-70,-95,-125,-Inf', 'Style', 'edit');
pard.binEdges.position = [8, 2.9];
pard.binEdges.Width = 2.1;
pard.binEdges.TooltipString = ['distance only: bin edges (nm), all decreasing or all increasing; bins appear left to right in this order, ' ...
    'each with its lower edge included. Default: >15, 15..0, 0..-10, ..., -95..-125, <-125.'];
pard.tLayout.object = struct('String', 'layout', 'Style', 'text');
pard.tLayout.position = [7, 1];
pard.layout.object = struct('String', {{'fixed', 'compact'}}, 'Style', 'popupmenu');
pard.layout.position = [7, 2];
pard.layout.Width = 0.8;
pard.layout.TooltipString = ['fixed: the fitted base of every site at (b x spacing, spacing), the same in every data set, so averages ' ...
    'of different data sets overlay exactly. compact: centred on the ROI centres, with gap nm between the averages (depends on the data).'];
pard.tSpacing.object = struct('String', 'spacing (nm)', 'Style', 'text');
pard.tSpacing.position = [7, 2.9];
pard.tSpacing.Width = 0.9;
pard.binSpacing.object = struct('String', '550', 'Style', 'edit');
pard.binSpacing.position = [7, 3.8];
pard.binSpacing.Width = 0.5;
pard.binSpacing.TooltipString = 'fixed layout: distance (nm) between neighbouring bin centres. Use the same value for all data sets you want to combine.';
pard.tGap.object = struct('String', 'gap (nm)', 'Style', 'text');
pard.tGap.position = [6, 1];
pard.tGap.Width = 1;
pard.binGap.object = struct('String', '250', 'Style', 'edit');
pard.binGap.position = [6, 2];
pard.binGap.Width = 0.5;
pard.binGap.TooltipString = ['compact layout only: empty space (nm) between neighbouring averages, edge to edge (edge: 95 % of the localizations). ' ...
    'Localizations further than half a bin spacing from their bin centre are left out, so bins never mix.'];

pard.tName.object = struct('String', 'name', 'Style', 'text');
pard.tName.position = [5, 1];
pard.name.object = struct('String', 'aligned_avg', 'Style', 'edit');
pard.name.position = [5, 2];
pard.writeRegister.object = struct('String', 'also write SMAP Register fields', 'Style', 'checkbox', 'Value', 0);
pard.writeRegister.position = [5, 3];
pard.writeRegister.Width = 2;
pard.writeRegister.TooltipString = ['After sorting, run all sites with Review only + Register: adds xnmaligned_/ynmaligned_/rank_ fields ' ...
    '(for mkMovieByFilter etc.). Not needed for the averages.'];

pard.plugininfo.name = 'LocMoFit_fitSortAverage';
pard.plugininfo.type = 'ROI_Analyze';
pard.plugininfo.description = ['Fit all sites with LocMoFit, sort by a fit parameter, align each site on its fit (model 1) ' ...
    'and average the sorted sites in N bins, from the first to the last usable site.'];
end
