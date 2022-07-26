%% Settings
rootPath = 'C:\Users\ries\git\SMAP\LocMoFit\rstGeneration\';
save2 = 'C:\Users\ries\git\SMAP\LocMoFit\docs\images\models\';
%% basic settings
par_general.rep = 1;
par_general.n = [1 1];
par_general.EMon = 1;
par_general.use_psf = 0;
par_general.saveTo = 'X:\users\Yu-Le\LocMoFit_t2\NPC3D\dualEllipseSim\simulation\';

% standard condition
par_std.LE = 0.7;
par_std.RB = 4;
par_std.linkageError = 3;
par_std.lifetime = 1;
par_std.frame = 1000000;
par_std.ellipticity = 0.1;

% brightness
par_mode{1}.name = 'storm';
par_mode{1}.model = 'Dye';
par_mode{1}.photon = 5000;
par_mode{1}.BG_pixel = 100;

%% Prepare common parameters
fn_mode = fieldnames(par_mode{1});

% Settings for the simulations
p = h_SimulateSites.getGuiParameters;

h_guiFile = g.children.guiFile;

p.labeling_efficiency = par_std.LE;
p.blinks = par_std.RB;
p.linkageerrorfree = par_std.linkageError;
p.lifetime = par_std.lifetime;
p.maxframes = par_std.frame;
p.numberofsites = par_general.n;
p.EMon = par_general.EMon;
p.use_psf = par_general.use_psf;
% p.psf_file = par_general.psf_file;

%% Get all the listed models
SMAP
locs.xnm = 0; locs.ynm = 0; locs.znm = 0; locs.locprecnm = 5; locs.locprecznm = 10; locs.layer = 1;
% SE
se = g.locData.SE;
eval_ = se.processors.eval;
h_SimulateSites = g.children.guiSites.children.Segment.children.SimulateSites;

roiSize = 200;
model2Show = modelList;
for k = 1:length(model2Show)
    eval_.addmodule('LocMoFitGUI');
    idxFitter = find(strcmp('LocMoFitGUI',fieldnames(eval_.children)));
    fitterGUI = eval_.processors{idxFitter};

    fitterGUI.load_callback([],[],[rootPath 'settings\modelPreview\' oneModel '_LocMoFit.mat'])
    fitterGUI.fitter.linkedGUI = fitterGUI;
    fitterGUI.updateGUI_fromLocMoFitObj;
    dimension = fitterGUI.fitter.model{1}.dimension;
    switch dimension
        case 2
            fig_pan = figure;
            pan = panel(fig_pan);
            pan.pack('h',2);
            ax = fitterGUI.fitter.plot(locs);
            fig_parent = ax.Parent;
            pan(1).select(ax);
            delete(findobj(ax,'type','line'));
            close(fig_parent)
        case 3
            fig_pan = figure;
            pan = panel(fig_pan);
            pan.pack('h',2);
            pan(1).pack('v',2);
            pan(2).pack('v',2);
    end
    
    
    inp = fitterGUI.getGuiParameters;
    inp.noFit = 0;
    fitterGUI.setGuiParameters(inp);

    h_SimulateSites.prepare_LocMoFit(fitterGUI.fitter);
    h_SimulateSites.setPar('useDepth', 0);
    h_SimulateSites.setPar('finalROISize','250');

    %% Perform simulations
    % settings are determined according to the k^th mode, the j^th varying
    % factor, and its l^th value/condition
    % For every brightness (mode)
    % For every varing factors (one of BG, LE, RB)
    % Reset all parameters

    % Set up brightness related parameters
    current_p = p;
    current_p.photons = par_mode{1}.photon;
    current_p.background = par_mode{1}.BG_pixel;
    mod_options = current_p.model.String;
    idx = find(strcmp(mod_options, par_mode{1}.model));
    current_p.model.Value = idx;
    current_p.model.selection = mod_options{idx};

    h_SimulateSites.setGuiParameters(current_p);
    h_SimulateSites.setPar('modelType', 'Point');
    for s = 1:par_general.rep
        rng(20220725)
        h_SimulateSites.processgo;
    end
    % clear files
    %             h_guiFile.call_menu_callback('clear')

    %% Make plot
    % set ROI
    h = g.getPar('sr_roihandle');
    pos = [500 500]./1000;
    rot = 90;
    linewidth_roi = 200;
    linelength_roi = 200;
    prepos = [0 -linelength_roi/2;0 linelength_roi/2]./1000;
    [prepos(:,1),prepos(:,2)] = rotcoord(prepos(:,1),prepos(:,2),deg2rad(rot));
    pos = prepos+pos;
    g.setPar('linewidth_roi',linewidth_roi);
    setPosition(h,pos)

    % Viewer3DV01
    module = g.children.guiAnalysis.children.sr3D.processors.Viewer3DV01;
    guiPar = module.getGuiParameters;
    guiPar.zmax = linewidth_roi/2;
    guiPar.zmin = -linewidth_roi/2;
    guiPar.pixrecset = [1 1];
    module.setGuiParameters(guiPar);

    module.guihandles.theta.String = '90';
    module.run(module.getAllParameters);

    % set layer
    L1 = g.children.guiRender.children.Layer1;
    guiPar = L1.getGuiParameters;
    guiPar.colorfield_min = 0;
    guiPar.colorfield_max = 1;
    guiPar.imax_min = -3.5;
    guiPar.lut.Value = 1;
    guiPar.render_colormode.Value = 1;
    L1.setGuiParameters(guiPar)

    guiFormat = g.getPar('guiFormat');
    guiFormat.guihandles.pixrec.String = '1';
    g.setPar('sr_colorbarthickness',0)
    g.setPar('sr_plotscalebar',0)

    L1.updateLayerField
    notify(g.P,'sr_render')

    switch dimension
        case 2
            % top data
            guiPar_L1_default = L1.getGuiParameters;
            oneView = copy(g.getPar('figurehandle_3Dviewer'));
            pan(2).select(oneView.Children);
            close(oneView);

            pan.margin = [0 0 0 5];
            pan.de.margin = [1 0 0 0];
            set(pan.de.axis, 'visible', 'on','xtick',[],'ytick',[],'xlabel',[],'ylabel',[])
            axis(pan.de.axis,'image')
            title(pan(1).axis, 'model');
            title(pan(2).axis, 'simulation');
        
            fig_pan.Position(3:4) = [500 265];

            sb = addScalebar(pan(2).axis,'bottom-right', roiSize./[15 15], 50);
            sb.LineWidth = 3;
        case 3
            % top data
            guiPar_L1_default = L1.getGuiParameters;
            oneView = copy(g.getPar('figurehandle_3Dviewer'));
            pan(2,1).select(oneView.Children);
            close(oneView);
        
            % side data
            guiPar.theta = 0;
            module.setGuiParameters(guiPar);
            module.run(module.getAllParameters);
            oneView = copy(g.getPar('figurehandle_3Dviewer'));
            pan(2,2).select(oneView.Children);
            close(oneView);
    end

    %% Save
    dateCreated = ['_c' char(datetime('now','TimeZone','local','Format','yyMMdd'))];
    mkdir([fPath_parent 'img\'])
    fName = [fPath_parent 'img\' inp.identifier '_singleSite' dateCreated '.pdf'];
    exportgraphics(fig_oneMT3D,fName,'ContentType','vector', 'Resolution', 600)

    eval_.removemodule('LocMoFitGUI');
end
%% Save
fName = [save2 oneModel '.png'];
exportgraphics(fig_pan,fName, 'Resolution', 200)
