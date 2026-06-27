classdef AverageSitesByDistance < interfaces.DialogProcessor & interfaces.SEProcessor
% Average sites grouped by distance bins (nm).
% Bin edges can be a single number (uniform bins of that width, starting
% from 0) or a comma-separated list of thresholds (custom edges, with
% -inf and +inf added automatically at both ends).
% Examples:
%   "20"        -> bins [0,20), [20,40), [40,60), ...
%   "10,20,50"  -> bins (-inf,10), [10,20), [20,50), [50,+inf)
    methods
        function obj = AverageSitesByDistance(varargin)
            obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters = {'se_viewer', 'se_siteroi', 'se_sitefov'};
            obj.showresults = true;
        end

        function out = run(obj, p)
            out = [];
            sites = obj.locData.SE.sites;
            if isempty(sites)
                obj.status('No sites loaded.'); return
            end
            roisize = ones(2,1) * obj.P.par.se_siteroi.content;

            % --- get per-site distance values ---
            distField = strtrim(p.distField);
            if isempty(distField)
                errordlg('Please specify a distance field.', 'AverageSitesByDistance');
                return
            end
            distances = getSiteDistances(sites, distField);
            validMask = ~isnan(distances);

            % apply distance range filter
            minD = str2double(num2str(p.minDist));
            maxD = str2double(num2str(p.maxDist));
            if isnan(minD), minD = 0;   end
            if isnan(maxD), maxD = Inf; end
            validMask = validMask & distances >= minD & distances <= maxD;

            if sum(validMask) == 0
                errordlg(['No valid values found in field: ' distField], 'AverageSitesByDistance');
                return
            end

            % --- parse bin edges ---
            binEdges = parseBinSpec(p.binSpec, distances(validMask));
            if isempty(binEdges) || length(binEdges) < 2
                errordlg('Could not parse bin specification.', 'AverageSitesByDistance');
                return
            end
            numBins = length(binEdges) - 1;

            % --- pre-count sites per bin (for histogram) ---
            binCounts = zeros(1, numBins);
            binLabels = cell(1, numBins);
            for l = 1:numBins
                lo = binEdges(l); hi = binEdges(l+1);
                if l < numBins
                    binCounts(l) = sum(distances >= lo & distances < hi);
                else
                    binCounts(l) = sum(distances >= lo & distances <= hi);
                end
                binLabels{l} = sprintf('[%s,%s)', fmt(lo), fmt(hi));
                obj.status(sprintf('Bin %d: [%s, %s) nm  -> %d sites', ...
                    l, fmt(lo), fmt(hi), binCounts(l)));
                drawnow
            end

            % --- set up averaging machinery (mirrors AverageSites_window) ---
            fdcal = figure(233);
            dcal = plugin('ROIManager', 'Evaluate', 'generalStatistics', fdcal, obj.P);
            dcal.attachLocData(obj.SE.locData);
            dcal.makeGui;

            xfirst = roisize(1);
            yfirst = roisize(2);
            x0 = 0; y0 = 0;
            locc = [];
            newfile = obj.locData.files.filenumberEnd;
            binLocs = {};   % per-bin: struct with xnm1,ynm1,xnm2,ynm2,label,nSites

            if p.inOneFile
                sitePerRow = max(1, ceil(numBins / max(1, p.rowSites)));
                flagPlotCol = 1;
                flagPlotRow = 1;
                newfile = obj.locData.files.filenumberEnd + 1;
                oneFilename = ['distBins_F' num2str(newfile)];
                if p.addfile
                    obj.locData.addfile(oneFilename);
                    initGuiAfterLoad(obj);
                    obj.SE.processors.preview.updateFilelist;
                end
            end

            for l = 1:numBins
                lo = binEdges(l); hi = binEdges(l+1);

                if l < numBins
                    inBin = distances >= lo & distances < hi;
                else
                    inBin = distances >= lo & distances <= hi;
                end
                binIdx = find(inBin);

                if isempty(binIdx)
                    obj.status(sprintf('Bin %d [%s,%s): 0 sites, skipping.', l, fmt(lo), fmt(hi)));
                    drawnow; continue
                end

                % --- apply per-bin site cap (for fast preview) ---
                maxS = round(p.maxSites);
                if maxS > 0 && length(binIdx) > maxS
                    binIdx = binIdx(randperm(length(binIdx), maxS));
                end

                if ~p.inOneFile
                    oneFilename = [fmt(lo) 'to' fmt(hi) 'nm'];
                    newfile = obj.locData.files.filenumberEnd + 1;
                    if p.addfile
                        obj.locData.addfile(oneFilename);
                        initGuiAfterLoad(obj);
                        obj.SE.processors.preview.updateFilelist;
                    end
                end

                if p.inOneFile
                    x0 = xfirst + (flagPlotCol-1) * roisize(1);
                    y0 = yfirst + (flagPlotRow-1) * roisize(2);
                    flagPlotCol = flagPlotCol + 1;
                    if flagPlotCol > sitePerRow
                        flagPlotRow = flagPlotRow + 1;
                        flagPlotCol = 1;
                    end
                end

                allIndLoc  = [];
                xnmrotAll  = [];
                ynmrotAll  = [];
                classAll   = [];
                ticc = tic;

                for ki = 1:length(binIdx)
                    k = binIdx(ki);
                    if p.sortselection.Value == 1 || sites(k).annotation.use
                        dcal.site = sites(k);
                        dcal.site.image = obj.locData.SE.plotsite(sites(k));
                        [locsSite, indloc] = dcal.getLocs({'xnmrot','ynmrot'}, ...
                            'size', roisize', 'grouping', 'ungrouped');
                        indloc = find(indloc);
                        if ~isempty(indloc)
                            allIndLoc(end+1:end+length(indloc), :)  = indloc;
                            xnmrotAll(end+1:end+length(indloc), :)  = locsSite.xnmrot;
                            ynmrotAll(end+1:end+length(indloc), :)  = locsSite.ynmrot;
                            classAll(end+1:end+length(indloc), :)   = sites(k).ID;
                        end
                    end
                    if toc(ticc) > 1
                        ticc = tic;
                        obj.status(sprintf('Bin %d/%d  site %d/%d', l, numBins, ki, length(binIdx)));
                        drawnow
                    end
                end

                if isempty(allIndLoc), continue, end

                locnew = obj.locData.loc;
                subsel = @(x) x(allIndLoc, :);
                locnew = structfun(subsel, locnew, 'UniformOutput', false);
                locnew.xnm       = xnmrotAll;
                locnew.ynm       = ynmrotAll;
                locnew.class     = classAll;
                locnew.origin    = allIndLoc;
                locnew.filenumber = repelem(newfile, length(allIndLoc))';

                locc = obj.locData.copy;
                locc.loc = locnew;
                locc.regroup;
                locc.filter;

                % --- store averaged locs for this bin (for subplot figure) ---
                layers_on = find(obj.getPar('sr_layerson'));
                lc1 = locc.getloc({'xnm','ynm'}, 'layer', 1, 'position', 'all', 'filenumber', newfile);
                lc2.xnm = []; lc2.ynm = [];
                if any(layers_on == 2)
                    lc2 = locc.getloc({'xnm','ynm'}, 'layer', 2, 'position', 'all', 'filenumber', newfile);
                end
                binLocs{end+1} = struct( ...
                    'xnm1', lc1.xnm, 'ynm1', lc1.ynm, ...
                    'xnm2', lc2.xnm, 'ynm2', lc2.ynm, ...
                    'label', binLabels{l}, 'nSites', length(binIdx));

                if p.addfile
                    locnew.xnm = locnew.xnm + x0;
                    locnew.ynm = locnew.ynm + y0;
                    fn = fieldnames(locnew);
                    for kf = 1:length(fn)
                        obj.locData.addloc(fn{kf}, locnew.(fn{kf}));
                    end
                    obj.locData.regroup;
                    obj.locData.filter;
                end
            end

            if isempty(locc)
                obj.status('No bins contained any sites.'); return
            end

            % --- Figure 1: histogram of sites per bin ---
            fhist = figure(901);
            clf(fhist);
            ax1 = axes(fhist);
            nonEmptyMask = binCounts > 0;
            bar(ax1, find(nonEmptyMask), binCounts(nonEmptyMask), 'FaceColor', [0.25 0.55 0.85]);
            set(ax1, 'XTick', 1:numBins, 'XTickLabel', binLabels, ...
                'XTickLabelRotation', 45);
            xlabel(ax1, 'Distance bin (nm)');
            ylabel(ax1, 'Number of sites');
            title(ax1, 'Sites per distance bin');
            grid(ax1, 'on');

            % --- Figure 2: averaged localization scatter, one subplot per bin ---
            if ~isempty(binLocs)
                nPlots  = length(binLocs);
                nCols   = min(nPlots, 5);
                nRows   = ceil(nPlots / nCols);
                hasTwo  = any(cellfun(@(b) ~isempty(b.xnm2), binLocs));

                favg = figure(902);
                clf(favg);
                set(favg, 'Name', 'Averaged sites per distance bin', 'NumberTitle', 'off');

                % global axis limits so all subplots are on the same scale
                allX = cellfun(@(b) [b.xnm1(:); b.xnm2(:)], binLocs, 'UniformOutput', false);
                allY = cellfun(@(b) [b.ynm1(:); b.ynm2(:)], binLocs, 'UniformOutput', false);
                allX = cat(1, allX{:});  allY = cat(1, allY{:});
                xLim = [min(allX) max(allX)];
                yLim = [min(allY) max(allY)];

                for bi = 1:nPlots
                    b  = binLocs{bi};
                    ax = subplot(nRows, nCols, bi, 'Parent', favg);
                    hold(ax, 'on');
                    plot(ax, b.xnm1, b.ynm1, '.', 'MarkerSize', 1, 'Color', [0.20 0.55 0.85]);
                    if hasTwo && ~isempty(b.xnm2)
                        plot(ax, b.xnm2, b.ynm2, '.', 'MarkerSize', 1, 'Color', [0.85 0.33 0.10]);
                    end
                    hold(ax, 'off');
                    axis(ax, 'equal');
                    xlim(ax, xLim);  ylim(ax, yLim);
                    box(ax, 'off');
                    set(ax, 'XTick', [], 'YTick', []);
                    title(ax, sprintf('%s nm  (n=%d)', b.label, b.nSites), ...
                        'FontSize', 7, 'Interpreter', 'none');
                end

                if hasTwo
                    sgtitle(favg, 'Averaged sites per distance bin  ●Ch1 (ring)  ●Ch2 (cap)');
                else
                    sgtitle(favg, 'Averaged sites per distance bin');
                end
            end
        end

        function pard = guidef(obj)
            pard = guidef(obj);
        end
    end
end

% -------------------------------------------------------------------------
% Button callbacks
% -------------------------------------------------------------------------

function selectField_callback(~, ~, obj)
% Browse site fields to fill the distance field edit box.
% If the selected field is a numeric array, shows an extra dialog
% listing each element together with its name (if a parallel .name
% field exists – as in allParsArg.value / allParsArg.name).
site = obj.SE.currentsite;
if isempty(site), return, end

str = browsefields(site, '');
if isempty(str) || strcmp(str, 'abortbutton'), return, end

% If the returned field is a numeric array we need an index
try
    val = eval(['site.' str]);
    if isnumeric(val) && numel(val) > 1
        % Try to find a companion .name field (e.g. allParsArg.name)
        parts = strsplit(str, '.');
        nameFieldParts      = parts;
        nameFieldParts{end} = 'name';
        nameFieldStr        = strjoin(nameFieldParts, '.');
        typeFieldParts      = parts;
        typeFieldParts{end} = 'type';
        typeFieldStr        = strjoin(typeFieldParts, '.');
        try
            nameVals = eval(['site.' nameFieldStr]);
            try, typeVals = eval(['site.' typeFieldStr]); catch, typeVals = {}; end
            listItems = cell(numel(val), 1);
            for k = 1:numel(val)
                nm = nameVals{k};
                if ~isempty(typeVals)
                    nm = [typeVals{k} '.' nm]; %#ok<AGROW>
                end
                listItems{k} = sprintf('%d:  %-25s = %.4g', k, nm, val(k));
            end
        catch
            listItems = arrayfun(@(k) sprintf('%d:  %.4g', k, val(k)), ...
                                 1:numel(val), 'UniformOutput', false)';
        end
        idx = listdlg('ListString', listItems, ...
                      'SelectionMode', 'single', ...
                      'Name',          'Select distance parameter', ...
                      'PromptString',  'Pick the parameter that encodes distance:', ...
                      'ListSize',      [340 300]);
        if isempty(idx), return, end
        str = [str '(' num2str(idx) ')'];
    end
catch
    % Cannot evaluate – use the string as-is
end

obj.guihandles.distField.String = str;
end

function previewBins_callback(~, ~, obj)
% Show bin counts without running the full averaging
sites = obj.locData.SE.sites;
if isempty(sites), disp('No sites'); return, end
p = obj.getAllParameters;
distField = strtrim(p.distField);
if isempty(distField), disp('No field specified'); return, end
distances = getSiteDistances(sites, distField);
validMask = ~isnan(distances);
if sum(validMask) == 0
    disp(['No valid values in field: ' distField]); return
end
binEdges = parseBinSpec(p.binSpec, distances(validMask));
if isempty(binEdges) || length(binEdges) < 2
    disp('Could not parse bin spec'); return
end
numBins = length(binEdges) - 1;
msg = sprintf('Distance field: %s\n\n', distField);
for l = 1:numBins
    lo = binEdges(l); hi = binEdges(l+1);
    if l < numBins
        n = sum(distances >= lo & distances < hi);
    else
        n = sum(distances >= lo & distances <= hi);
    end
    msg = [msg sprintf('  Bin %d: [%s, %s) nm  ->  %d sites\n', ...
        l, fmt(lo), fmt(hi), n)]; %#ok<AGROW>
end
msg = [msg sprintf('\nTotal valid sites: %d / %d', sum(validMask), length(sites))];
msgbox(msg, 'Bin preview', 'help');
end

% -------------------------------------------------------------------------
% Local helpers
% -------------------------------------------------------------------------

function distances = getSiteDistances(sites, distField)
% Evaluate distField for every site.
%
% Two syntaxes are supported:
%
%  1. Field path  (default):
%       evaluation.LocMoFitGUI.allParsArg.value(14)
%     Evaluated as  s.<distField>  where s = sites(k).
%
%  2. Formula  (prefix with @):
%       @sqrt(v(1)^2 + (v(2)-v(14))^2)
%     The expression after @ is evaluated with these variables in scope:
%       s  = sites(k)                           (full site struct)
%       ap = s.evaluation.LocMoFitGUI.allParsArg  (allParsArg struct, if present)
%       v  = ap.value                           (parameter value vector)
%     Use this when the distance must be computed from multiple parameters,
%     e.g. the Euclidean distance between a ring and a cap model.
isFormula = ischar(distField) && ~isempty(distField) && distField(1) == '@';
if isFormula
    expr = distField(2:end);   % strip the leading '@'
end

% debug: print what we received
fprintf('getSiteDistances: ischar=%d  isFormula=%d  distField="%s"\n', ...
    ischar(distField), isFormula, num2str(distField));

distances = NaN(1, length(sites));
firstErr = '';
for k = 1:length(sites)
    try
        s = sites(k); %#ok<NASGU>
        if isFormula
            % make ap and v available if allParsArg exists
            try
                ap = s.evaluation.LocMoFitGUI.allParsArg; %#ok<NASGU>
                v  = ap.value; %#ok<NASGU>
            catch
            end
            val = eval(expr);
        else
            val = eval(['s.' distField]);
            if isstruct(val)
                val = val.value;
            end
            val = val(1);
        end
        distances(k) = val;
    catch err
        if isempty(firstErr)
            firstErr = err.message;
        end
    end
end
if ~isempty(firstErr)
    fprintf('getSiteDistances: first error was: %s\n', firstErr);
end
end

function edges = parseBinSpec(specStr, validDists)
% Returns bin edges as a row vector.
% Single number  -> uniform bins starting from 0.
% Multiple numbers -> custom edges with -inf and +inf appended.
%
% Note: SMAP's getAllParameters() converts single-number edit-box strings
% to doubles automatically, so specStr may arrive as a numeric scalar.
if isnumeric(specStr)
    nums = specStr;
else
    specStr = strtrim(specStr);
    parts = strsplit(specStr, ',');
    nums = cellfun(@str2double, strtrim(parts));
    nums = nums(~isnan(nums));
end

if isempty(nums)
    edges = [];
    return
end

if isscalar(nums)
    bw = nums;
    if bw <= 0, edges = []; return, end
    hi = ceil(max(validDists) / bw) * bw;
    edges = 0 : bw : hi;
    if isempty(edges) || edges(end) < max(validDists)
        edges(end+1) = edges(end) + bw;
    end
else
    % Custom thresholds: always start from 0, end at +inf
    edges = [0, sort(nums(nums > 0)), inf];
end
end

function s = fmt(v)
% Format a bin edge value for display / filename
if isinf(v) && v < 0
    s = '-inf';
elseif isinf(v)
    s = 'inf';
else
    s = num2str(v);
end
end

% -------------------------------------------------------------------------
% GUI definition
% -------------------------------------------------------------------------

function pard = guidef(obj)

pard.tDist.object  = struct('String', 'Distance field', 'Style', 'text');
pard.tDist.position = [1, 1];
pard.tDist.Width   = 1;

pard.distField.object  = struct('String', '', 'Style', 'edit');
pard.distField.position = [1, 2];
pard.distField.Width   = 2;
pard.distField.TooltipString = [...
    'Path to a per-site scalar, e.g.  evaluation.LocMoFitGUI.allParsArg.value(14)' char(10) ...
    'For computed distances prefix with @, e.g.:' char(10) ...
    '  @sqrt(v(1)^2+(v(2)-v(14))^2)' char(10) ...
    'Inside @ expressions: s=site, ap=allParsArg, v=allParsArg.value'];

pard.selectField.object = struct('String', 'select', 'Style', 'pushbutton', ...
    'Callback', {{@selectField_callback, obj}});
pard.selectField.position = [1, 4];
pard.selectField.Width = 0.7;

pard.tBin.object  = struct('String', 'Bin width or edges (nm)', 'Style', 'text');
pard.tBin.position = [2, 1];
pard.tBin.Width   = 1.5;

pard.binSpec.object  = struct('String', '20', 'Style', 'edit');
pard.binSpec.position = [2, 2.5];
pard.binSpec.Width   = 1.5;
pard.binSpec.TooltipString = ['Single number: uniform bin width starting from 0.  e.g. "20" -> [0,20), [20,40), ...' char(10) ...
    'Comma-separated: custom edges from 0.  e.g. "10,70,200" -> [0,10), [10,70), [70,200), [200,+inf)'];

pard.previewBins.object = struct('String', 'preview', 'Style', 'pushbutton', ...
    'Callback', {{@previewBins_callback, obj}});
pard.previewBins.position = [2, 4];
pard.previewBins.Width = 0.7;

pard.tSort.object  = struct('String', 'Sites to use', 'Style', 'text');
pard.tSort.position = [3, 1];
pard.tSort.Width   = 1;

pard.sortselection.object = struct('String', {{'all', 'use'}}, 'Style', 'popupmenu');
pard.sortselection.position = [3, 2];
pard.sortselection.Width   = 1;

pard.tName.object  = struct('String', 'Name', 'Style', 'text');
pard.tName.position = [4, 1];
pard.tName.Width   = 1;

pard.name.object  = struct('String', 'average', 'Style', 'edit');
pard.name.position = [4, 2];
pard.name.Width   = 1;

pard.addfile.object  = struct('String', 'add averages as new datasets', 'Style', 'checkbox', 'Value', 1);
pard.addfile.position = [4, 3];
pard.addfile.Width   = 2;

pard.tRow.object  = struct('String', 'Bins per row', 'Style', 'text');
pard.tRow.position = [5, 1];
pard.tRow.Width   = 1;

pard.rowSites.object  = struct('String', '4', 'Style', 'edit');
pard.rowSites.position = [5, 2];
pard.rowSites.Width   = 0.5;

pard.inOneFile.object  = struct('String', 'all in one file', 'Style', 'checkbox', 'Value', 1);
pard.inOneFile.position = [5, 3];
pard.inOneFile.Width   = 1.5;

pard.tMaxSites.object  = struct('String', 'Max sites/bin (0=all)', 'Style', 'text');
pard.tMaxSites.position = [6, 1];
pard.tMaxSites.Width   = 1.5;

pard.maxSites.object  = struct('String', '0', 'Style', 'edit');
pard.maxSites.position = [6, 2.5];
pard.maxSites.Width   = 0.5;

pard.tDistRange.object  = struct('String', 'Distance range (nm)', 'Style', 'text');
pard.tDistRange.position = [7, 1];
pard.tDistRange.Width   = 1.5;

pard.minDist.object  = struct('String', '0', 'Style', 'edit');
pard.minDist.position = [7, 2.5];
pard.minDist.Width   = 0.5;
pard.minDist.TooltipString = 'Minimum distance (nm). Sites below this are excluded.';

pard.tDistTo.object  = struct('String', 'to', 'Style', 'text');
pard.tDistTo.position = [7, 3];
pard.tDistTo.Width   = 0.3;

pard.maxDist.object  = struct('String', 'Inf', 'Style', 'edit');
pard.maxDist.position = [7, 3.3];
pard.maxDist.Width   = 0.5;
pard.maxDist.TooltipString = 'Maximum distance (nm). Sites above this are excluded. Leave as Inf for no upper limit.';

pard.plugininfo.name = 'AverageSitesByDistance';
pard.plugininfo.type = 'ROI_Analyze';
pard.plugininfo.description = ...
    'Average sites grouped by a per-site distance parameter. Bin by uniform width or custom nm thresholds.';

end
