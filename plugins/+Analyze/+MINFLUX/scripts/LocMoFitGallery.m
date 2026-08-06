% Gallery of the LocMoFit fits: for every site that has a fit image, the site
% image and the site image overlaid with the fitted model are plotted next to
% each other on A4 pages.
% The fit image is written by the LocMoFitGUI evaluation plugin. It is rendered
% on the same grid as the site image, so the two can be overlaid pixel by pixel.

% parameters
p.evalname = '';            % name of the LocMoFit evaluation, empty: take the first one found
p.sitesPerRow = 3;          % number of site/overlay pairs next to each other
p.rowsPerPage = 5;          % number of rows on one A4 page
p.lutModel = {'green'};     % lut for the fitted model, one per layer, cycled
p.modelWeight = 1;          % how strong the model is added on top of the site image
p.modelSaturation = 0.7;    % fraction of the maximum that is shown white, smaller: more contrast (and clipping)
p.modelGamma = .8;         % <1 brightens the weak parts of the model
p.dataWeight = 0.7;         % the site image is dimmed by this factor in the overlay
p.renderMissing = true;     % render the site image if it has not been rendered yet
p.showID = true;            % write the site ID above each pair
p.pdffile = '';             % if set, all pages are appended into this pdf

%%
SE = g.locData.SE;
sites = SE.sites;

if isempty(p.evalname)
    for k = 1:length(sites)
        fn = fieldnames(sites(k).evaluation);
        % only an evaluation that actually wrote a fit image is of use here
        fn = fn(startsWith(fn, 'LocMoFit'));
        for l = 1:length(fn)
            if isfield(sites(k).evaluation.(fn{l}), 'image')&&~isempty(sites(k).evaluation.(fn{l}).image)
                p.evalname = fn{l};
                break
            end
        end
        if ~isempty(p.evalname)
            break
        end
    end
    if isempty(p.evalname)
        error('No LocMoFit evaluation found. Please set p.evalname.')
    end
    disp(['Using the evaluation ' p.evalname '.'])
end

% collect the sites that have a fit image
indUse = [];
for k = 1:length(sites)
    if ~isfield(sites(k).evaluation, p.evalname)
        continue
    end
    % sites evaluated before the fit image was introduced have no image field
    if ~isfield(sites(k).evaluation.(p.evalname), 'image')
        continue
    end
    imfit = sites(k).evaluation.(p.evalname).image;
    if isstruct(imfit)&&~isempty(imfit.image)
        indUse(end+1) = k; %#ok<SAGROW>
    end
end
if isempty(indUse)
    error(['None of the sites has a fit image. Re-run the evaluation with ' p.evalname '.'])
end
disp([num2str(length(indUse)) ' of ' num2str(length(sites)) ' sites have a fit image.'])

perPage = p.sitesPerRow*p.rowsPerPage;
numPages = ceil(length(indUse)/perPage);
if ~isempty(p.pdffile)&&exist(p.pdffile,'file')
    delete(p.pdffile)
end

for pg = 1:numPages
    f = figure('Units','centimeters','Position',[2 2 21 29.7],'Color','w',...
        'PaperUnits','centimeters','PaperSize',[21 29.7],'PaperPosition',[0 0 21 29.7],...
        'Name',['LocMoFit gallery ' num2str(pg) '/' num2str(numPages)]);
    % two tiles per site: the site and the site with the model on top
    t = tiledlayout(f, p.rowsPerPage, 2*p.sitesPerRow,'TileSpacing','tight','Padding','compact');
    title(t, [p.evalname ' (' num2str(pg) '/' num2str(numPages) ')'],'Interpreter','none')

    indPage = indUse((pg-1)*perPage+1:min(pg*perPage, length(indUse)));
    for k = 1:length(indPage)
        site = sites(indPage(k));
        if (~isstruct(site.image)||isempty(site.image))&&p.renderMissing
            SE.plotsite(site);
        end
        if ~isstruct(site.image)||isempty(site.image)
            warning(['Site ' num2str(site.ID) ' has no site image, skipped.'])
            continue
        end
        imSite = im2double(site.image.image);
        imModel = modelRGB(site.evaluation.(p.evalname).image, p.lutModel, p.modelSaturation, p.modelGamma);
        if ~isequal(size(imModel,1,2), size(imSite,1,2))
            % should not happen, the fit is rendered on the grid of the site image
            warning(['Site ' num2str(site.ID) ': the fit and the site image differ in size, the fit is resized.'])
            imModel = imresize(imModel, size(imSite,1,2));
        end
        imOverlay = min(p.dataWeight*imSite + p.modelWeight*imModel, 1);

        ax1 = nexttile(t);
        showimage(ax1, imSite);
        if p.showID
            title(ax1, ['site ' num2str(site.ID)],'FontSize',7,'FontWeight','normal')
        end
        ax2 = nexttile(t);
        showimage(ax2, imOverlay);
        if p.showID
            title(ax2, 'fit','FontSize',7,'FontWeight','normal')
        end
    end

    if ~isempty(p.pdffile)
        exportgraphics(f, p.pdffile, 'Append', true)
    end
end
if ~isempty(p.pdffile)
    disp(['Gallery written to ' p.pdffile '.'])
end

function im = modelRGB(imfit, allLut, saturation, gammaVal)
% the fit image holds intensities only, one plane per layer: stretch the
% contrast and apply a lut here
im = 0;
for k = 1:size(imfit.image,3)
    val = double(imfit.image(:,:,k))/255;
    val = min(val./saturation, 1).^gammaVal;    % clip the top and lift the weak parts
    lut = mymakelut(allLut{mod(k-1,length(allLut))+1});
    im = im + ind2rgb(uint8(val*255), lut);
end
im = min(im, 1);
end

function showimage(ax, im)
image(ax, im)
axis(ax, 'image')      % keep the aspect ratio
ax.XTick = [];
ax.YTick = [];
ax.Box = 'on';
end
