classdef SiteAutoDetect<interfaces.DialogProcessor&interfaces.SEProcessor
% Automatically detects clathrin-coated pit (CCP) sites from SMLM
% localizations and adds them to the SMAP Site Explorer.
%
% Algorithm: density map → Gaussian blur → local maxima → min-distance
% suppression → centroid refinement → optional locs-count filter → SE.
    methods
        function obj=SiteAutoDetect(varargin)
            obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'se_cellfov','se_siteroi'};
        end

        function out=run(obj,p)
            params = obj.getParams();
            if isempty(params), out=[]; return; end
            obj.doDetect(params, false);
            out=[];
        end

        function initGui(obj)
            obj.guihandles.preview_button.Callback = @(~,~) obj.runPreview();
            obj.guihandles.detect_button.Callback  = @(~,~) obj.runDetect();
        end

        function runPreview(obj)
            params = obj.getParams();
            if isempty(params), return; end
            obj.doDetect(params, true);
        end

        function runDetect(obj)
            params = obj.getParams();
            if isempty(params), return; end
            obj.doDetect(params, false);
        end

        function pard=guidef(obj)
            pard=guidef_SiteAutoDetect(obj);
        end
    end

    methods(Access=private)
        function params = getParams(obj)
            params = [];
            try
                params.threshold    = str2double(obj.guihandles.threshold_edit.String);
                params.min_dist_nm  = str2double(obj.guihandles.min_dist_edit.String);
                params.blur_sigma   = str2double(obj.guihandles.blur_sigma_edit.String);
                params.min_locs     = str2double(obj.guihandles.min_locs_edit.String);
                params.locs_bins    = str2double(obj.guihandles.locs_bins_edit.String);
                params.clear_sites  = obj.guihandles.clear_sites.Value;
                strs = obj.guihandles.filter_mode.String;
                params.filter_mode  = strs{obj.guihandles.filter_mode.Value}; % 'Exclude' or 'Mark'
            catch
                warndlg('Could not read GUI parameters.'); return;
            end
            if any(isnan([params.threshold params.min_dist_nm params.blur_sigma params.min_locs]))
                warndlg('Invalid parameter values.'); params=[]; return;
            end
            % Read SE parameters
            try
                params.cell_fov  = obj.getPar('se_cellfov');
                params.site_roi  = obj.getPar('se_siteroi');
            catch
                params.cell_fov  = 5000;
                params.site_roi  = 300;
            end
        end

        function doDetect(obj, params, preview_only)
            threshold   = params.threshold;
            min_dist_nm = params.min_dist_nm;
            blur_sigma  = params.blur_sigma;
            min_locs    = params.min_locs;
            site_roi    = params.site_roi;
            cell_fov    = params.cell_fov;
            pixel_nm    = 10;

            % ── 1. Get all localizations ──────────────────────────────────
            try
                locs = obj.locData.getloc({'xnm','ynm'}, 'removeFilter', 'filenumber');
            catch
                locs = obj.locData.getloc({'xnm','ynm'});
            end
            xnm = double(locs.xnm(:));
            ynm = double(locs.ynm(:));
            if isempty(xnm)
                warndlg('No localizations loaded.'); return;
            end

            % ── 2. Build density map ──────────────────────────────────────
            x_min = min(xnm);  y_min = min(ynm);
            W = ceil((max(xnm) - x_min) / pixel_nm) + 1;
            H = ceil((max(ynm) - y_min) / pixel_nm) + 1;
            xi = min(W, max(1, round((xnm - x_min) / pixel_nm) + 1));
            yi = min(H, max(1, round((ynm - y_min) / pixel_nm) + 1));
            img = accumarray([yi, xi], 1, [H W]);

            % ── 3. Blur ───────────────────────────────────────────────────
            sigma_px = blur_sigma / pixel_nm;
            if exist('imgaussfilt','file') || exist('imgaussfilt','builtin')
                img_blur = imgaussfilt(double(img), sigma_px);
            else
                h = fspecial('gaussian', 2*ceil(3*sigma_px)+1, sigma_px);
                img_blur = imfilter(double(img), h, 'replicate');
            end

            % 99th-percentile normalisation for better contrast in preview
            p99 = prctile(img_blur(img_blur>0), 99);
            if p99 == 0, p99 = max(img_blur(:)) + eps; end
            img_norm = img_blur / p99;
            img_norm(img_norm > 1) = 1;

            % ── 4. Find local maxima ──────────────────────────────────────
            thresh_abs = threshold;   % threshold is already on 0-1 scale of img_norm
            if exist('imregionalmax','file') || exist('imregionalmax','builtin')
                bw = imregionalmax(img_norm) & (img_norm > thresh_abs);
            else
                se_dil = strel('disk', max(1, round(sigma_px)));
                bw = (imdilate(img_norm, se_dil) == img_norm) & (img_norm > thresh_abs);
            end
            [rows, cols] = find(bw);
            if isempty(rows)
                msgbox('No sites detected. Try lowering the threshold.'); return;
            end

            % ── 5. To nm, sort, min-distance suppression ─────────────────
            x_peaks = (cols - 1) * pixel_nm + x_min;
            y_peaks = (rows - 1) * pixel_nm + y_min;
            intensities = img_norm(sub2ind(size(img_norm), rows, cols));
            [~, order]  = sort(intensities, 'descend');
            x_peaks = x_peaks(order);
            y_peaks = y_peaks(order);

            N = length(x_peaks);
            kept = true(N, 1);
            if N > 1000 && (exist('knnsearch','file') || exist('knnsearch','builtin'))
                coords = [x_peaks, y_peaks];
                for k = 1:N
                    if ~kept(k), continue; end
                    idx = knnsearch(coords, coords(k,:), 'K', min(N,200));
                    idx = idx(idx > k);
                    d   = sqrt(sum((coords(idx,:) - coords(k,:)).^2, 2));
                    kept(idx(d < min_dist_nm)) = false;
                end
            else
                for k = 1:N
                    if ~kept(k), continue; end
                    dists = sqrt((x_peaks - x_peaks(k)).^2 + (y_peaks - y_peaks(k)).^2);
                    dists(k) = Inf;
                    rm = find(kept & dists < min_dist_nm & (1:N)' > k);
                    kept(rm) = false;
                end
            end
            x_peaks = x_peaks(kept);
            y_peaks = y_peaks(kept);

            % ── 6. Centroid refinement (uses se_siteroi) ──────────────────
            filenumber = [];
            if ~isempty(obj.SE.files) && obj.SE.numberOfFiles > 0
                filenumber = obj.SE.files(1).ID;
            end
            for k = 1:length(x_peaks)
                try
                    if ~isempty(filenumber)
                        rl = obj.locData.getloc({'xnm','ynm'}, ...
                            'position', [x_peaks(k), y_peaks(k), site_roi], ...
                            'filenumber', filenumber);
                    else
                        rl = obj.locData.getloc({'xnm','ynm'}, ...
                            'position', [x_peaks(k), y_peaks(k), site_roi]);
                    end
                    if ~isempty(rl.xnm) && length(rl.xnm) >= 3
                        x_peaks(k) = mean(double(rl.xnm));
                        y_peaks(k) = mean(double(rl.ynm));
                    end
                catch
                end
            end

            % ── 7. Count locs per site & apply min_locs filter ───────────
            n_sites = length(x_peaks);
            locs_count = zeros(n_sites, 1);
            for k = 1:n_sites
                try
                    if ~isempty(filenumber)
                        rl = obj.locData.getloc({'xnm'}, ...
                            'position', [x_peaks(k), y_peaks(k), site_roi], ...
                            'filenumber', filenumber);
                    else
                        rl = obj.locData.getloc({'xnm'}, ...
                            'position', [x_peaks(k), y_peaks(k), site_roi]);
                    end
                    locs_count(k) = length(rl.xnm);
                catch
                    locs_count(k) = 0;
                end
            end

            passes = locs_count >= min_locs;

            % ── 8. Update GUI label ───────────────────────────────────────
            n_found = sum(passes);
            try
                obj.guihandles.nsites_label.String = sprintf('N sites found: %d  (excl. %d by locs filter)', ...
                    n_found, sum(~passes));
            catch
            end

            % ── 9. Preview ────────────────────────────────────────────────
            if preview_only
                cols_k = round((x_peaks - x_min) / pixel_nm) + 1;
                rows_k = round((y_peaks - y_min) / pixel_nm) + 1;

                figure('Name','SiteAutoDetect Preview','Color','k', ...
                    'Units','normalized','Position',[0.1 0.1 0.8 0.45]);

                % Left: density map
                subplot(1,2,1);
                imagesc(img_norm); colormap(gca, hot); axis equal tight; colorbar;
                hold on;
                % Passed sites: cyan circles
                plot(cols_k(passes),  rows_k(passes),  'co', 'MarkerSize',10,'LineWidth',2);
                % Filtered sites: yellow x
                plot(cols_k(~passes), rows_k(~passes), 'yx', 'MarkerSize', 8,'LineWidth',1.5);
                title(sprintf('threshold=%.3f  |  %d sites  (%d below locs filter)', ...
                    threshold, n_found, sum(~passes)), 'Color','w');
                legend({'pass','locs<min'},'TextColor','w','Color','none','EdgeColor','none');

                % Right: histogram of locs per site
                subplot(1,2,2);
                bins = max(1, round(params.locs_bins));
                if ~isempty(locs_count)
                    histogram(locs_count, bins, 'FaceColor',[0.2 0.6 1], 'EdgeColor','none');
                end
                hold on;
                yl = ylim;
                plot([min_locs min_locs], yl, 'r--', 'LineWidth', 2);
                xlabel('Locs per site'); ylabel('N ROIs');
                title('Locs distribution', 'Color','w');
                set(gca,'Color','k','XColor','w','YColor','w');
                legend({sprintf('histogram'), sprintf('min\\_locs=%.0f',min_locs)}, ...
                    'TextColor','w','Color','none','EdgeColor','none');

                set(gcf,'Color','k');
                fprintf('[SiteAutoDetect] Preview: %d/%d sites pass, locs filter=%d\n', ...
                    n_found, length(x_peaks), min_locs);
                return;
            end

            % ── 10. Add to SE ─────────────────────────────────────────────
            se = obj.SE;
            if isempty(se.files) || se.numberOfFiles == 0
                warndlg('No files in Site Explorer. Load data first.'); return;
            end
            filenumber = se.files(1).ID;

            % Clear existing sites if requested
            if params.clear_sites && se.numberOfSites > 0
                ids = [se.sites(:).ID];
                for k = 1:length(ids)
                    try se.removeSite(ids(k)); catch, end
                end
                if se.numberOfCells > 0
                    cids = [se.cells(:).ID];
                    for k = 1:length(cids)
                        try se.removeCell(cids(k)); catch, end
                    end
                end
            end

            % Build cell grid: one cell per grid tile that contains ≥1 site
            half = cell_fov / 2;
            x_fov_min = min(x_peaks) - half;
            y_fov_min = min(y_peaks) - half;
            % Grid indices for each site
            gx = floor((x_peaks - x_fov_min) / cell_fov);  % 0-based
            gy = floor((y_peaks - y_fov_min) / cell_fov);
            grid_ids = gx * 1e6 + gy;  % unique key per tile

            [unique_tiles, ~, tile_idx] = unique(grid_ids);
            cellID_map = containers.Map('KeyType','double','ValueType','int32');

            for t = 1:length(unique_tiles)
                tile = unique_tiles(t);
                gxi = mod(tile, 1e6);   % recover gx
                gyi = tile - gxi * 1e6; % recover gy  (note: tile=gx*1e6+gy)
                cx = x_fov_min + (gxi + 0.5) * cell_fov;
                cy = y_fov_min + (gyi + 0.5) * cell_fov;
                newcell = interfaces.SEsites;
                newcell.pos = [cx, cy, 0];
                newcell.ID  = 0;
                newcell.info.filenumber = filenumber;
                cid = se.addCell(newcell);
                cellID_map(double(tile)) = int32(cid);
            end

            % Add sites
            n_added = 0;
            for k = 1:length(x_peaks)
                use_flag = passes(k);
                if strcmp(params.filter_mode, 'Exclude') && ~use_flag
                    continue;
                end
                site = interfaces.SEsites;
                site.pos  = [x_peaks(k), y_peaks(k), 0];
                site.ID   = 0;
                site.info.filenumber = filenumber;
                site.info.cell       = double(cellID_map(grid_ids(k)));
                site.annotation.use  = use_flag;
                se.addSite(site);
                n_added = n_added + 1;
            end

            % Refresh SE display
            try
                se.processors.preview.updateCelllist;
                se.processors.preview.updateSitelist;
            catch
            end

            fprintf('[SiteAutoDetect] Added %d sites to SE (%d marked Use=false by locs filter)\n', ...
                n_added, sum(~passes & strcmp(params.filter_mode,'Mark')));
        end
    end
end


function pard = guidef_SiteAutoDetect(obj) %#ok<INUSD>

% ── Row 1: Threshold ─────────────────────────────────────────────────────
pard.t_thresh.object = struct('String','Threshold (0–1)','Style','text');
pard.t_thresh.position = [1,1];
pard.t_thresh.Width = 1.2;

pard.threshold_edit.object = struct('String','0.02','Style','edit');
pard.threshold_edit.position = [1,2.3];
pard.threshold_edit.Width = 0.6;

pard.t_thresh_hint.object = struct('String', ...
    'Fraction of 99th-percentile map value. Lower = more (faint) sites.', ...
    'Style','text');
pard.t_thresh_hint.position = [2,1];
pard.t_thresh_hint.Width = 3;

% ── Row 3: Min site distance ──────────────────────────────────────────────
pard.t_dist.object = struct('String','Min site distance (nm)','Style','text');
pard.t_dist.position = [3,1];
pard.t_dist.Width = 1.2;

pard.min_dist_edit.object = struct('String','150','Style','edit');
pard.min_dist_edit.position = [3,2.3];
pard.min_dist_edit.Width = 0.6;

% ── Row 4: Blur sigma ────────────────────────────────────────────────────
pard.t_blur.object = struct('String','Blur sigma (nm)','Style','text');
pard.t_blur.position = [4,1];
pard.t_blur.Width = 1.2;

pard.blur_sigma_edit.object = struct('String','40','Style','edit');
pard.blur_sigma_edit.position = [4,2.3];
pard.blur_sigma_edit.Width = 0.6;

% ── Row 5: Min locs per site ─────────────────────────────────────────────
pard.t_minlocs.object = struct('String','Min locs per site','Style','text');
pard.t_minlocs.position = [5,1];
pard.t_minlocs.Width = 1.2;

pard.min_locs_edit.object = struct('String','10','Style','edit');
pard.min_locs_edit.position = [5,2.3];
pard.min_locs_edit.Width = 0.6;

% ── Row 6: Filter mode + histogram bins ──────────────────────────────────
pard.t_filtermode.object = struct('String','Locs filter action','Style','text');
pard.t_filtermode.position = [6,1];
pard.t_filtermode.Width = 1.2;

pard.filter_mode.object = struct('String',{{'Exclude','Mark (Use=false)'}},'Style','popupmenu');
pard.filter_mode.position = [6,2.3];
pard.filter_mode.Width = 1.2;

pard.t_bins.object = struct('String','Histogram bins','Style','text');
pard.t_bins.position = [7,1];
pard.t_bins.Width = 1.2;

pard.locs_bins_edit.object = struct('String','30','Style','edit');
pard.locs_bins_edit.position = [7,2.3];
pard.locs_bins_edit.Width = 0.6;

% ── Row 8: Clear checkbox ────────────────────────────────────────────────
pard.clear_sites.object = struct('String','Clear existing sites first','Style','checkbox','Value',0);
pard.clear_sites.position = [8,1];
pard.clear_sites.Width = 2;

% ── Row 9: Buttons ───────────────────────────────────────────────────────
pard.preview_button.object = struct('String','Preview','Style','pushbutton');
pard.preview_button.position = [9,1];
pard.preview_button.Width = 1;

pard.detect_button.object = struct('String','Detect Sites','Style','pushbutton');
pard.detect_button.position = [9,2];
pard.detect_button.Width = 1;

% ── Row 10: Status label ─────────────────────────────────────────────────
pard.nsites_label.object = struct('String','N sites found: ---','Style','text');
pard.nsites_label.position = [10,1];
pard.nsites_label.Width = 3;

pard.plugininfo.name        = 'Auto-detect CCP Sites';
pard.plugininfo.description = ['Detects CCP sites from SMLM localizations via density map + local maxima. ' ...
    'Threshold is a fraction of the 99th-percentile map intensity. ' ...
    'Works for both 2D and 3D data.'];
pard.plugininfo.type        = 'ProcessorPlugin';
end
