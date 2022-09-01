%% Prepare common parameters
% Settings for the simulations
SMAP
p = [];

%% basic settings
p.general.rep = 1;
p.general.saveTo = 'C:\Users\ries\git\SMAP\LocMoFit\docs\images\models\';
p.general.rootPath = 'C:\Users\ries\git\SMAP\LocMoFit\img_model\';
p.general.crop = [0 0]; % from both sides
p.general.thichness = par_general.roiSize; % cross-section
p.general.view = {'xy'}; % top-view (xy) or side-view (xz, yz)

h_SimulateSites = g.children.guiSites.children.Segment.children.SimulateSites;
p.SimulateSites = h_SimulateSites.getGuiParameters;

p.SimulateSites.labeling_efficiency = 0.7;
p.SimulateSites.blinks = 4;
p.SimulateSites.linkageerrorfree = 5;
p.SimulateSites.lifetime = 1;
p.SimulateSites.maxframes = 1000000;
p.SimulateSites.numberofsites = [1 1];
p.SimulateSites.EMon = 1;
p.SimulateSites.use_psf = 0;
p.SimulateSites.model.selection = 'dye';
p.SimulateSites.photons = 5000;
p.SimulateSites.background = 100;
% p.psf_file = par_general.psf_file;

h_settings = g.children.guiSites.children.settings;
p.settings = h_settings.getGuiParameters;
p.settings.se_sitefov = 200;
p.settings.se_siteroi = 200;

h_guiFile = g.children.guiFile;
%% Get all the listed models
locs.xnm = 0; locs.ynm = 0; locs.znm = 0; locs.locprecnm = 5; locs.locprecznm = 10; locs.layer = 1;
% SE

model_imgDefined = img_creation_settings();
for k = 1:length(model_imgDefined)
    oneModel = model_imgDefined{k};
    p_model = img_creation_settings(oneModel);

    fitterGUI = setUpLocMoFit(g, oneModel);
    modPar(fitterGUI, h_SimulateSites, h_settings, p, p_model);

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
    g.locData.setPar('layer1_selectedField', {'locprecnm', 0,10, 1, 1})

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

            pan.margin = [0 0 0 7];
            pan.de.margin = [1 0 0 0];
            pan.fontsize = 16;
            set(pan.de.axis, 'visible', 'on','xtick',[],'ytick',[],'xlabel',[],'ylabel',[])
            axis(pan.de.axis,'image')
            title(pan(1).axis, 'model');
            title(pan(2).axis, 'simulation');
        
            fig_pan.Position(3:4) = [500 274];

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
fName = [par_general.saveTo oneModel '.png'];
exportgraphics(fig_pan,fName, 'Resolution', 200)

%% Functions
function fitterGUI = setUpLocMoFit(g, oneModel)
    se = g.locData.SE;
    eval_ = se.processors.eval;
    eval_.addmodule('LocMoFitGUI');
    idxFitter = find(strcmp('LocMoFitGUI',fieldnames(eval_.children)));
    fitterGUI = eval_.processors{idxFitter};
    fitterGUI.fitter.addModel(functionModel(modelList(oneModel)))
    fitterGUI.fitter.linkedGUI = fitterGUI;
    fitterGUI.updateGUI_fromLocMoFitObj;
end

function modPar(fitterGUI, h_SimulateSites, h_settings, p, p_model)
% 20220818: need to work on this part
    fn = fieldnames(p_model);
    for k = 1:length(fn)
        p_SMAP = p.(fn{k}).SMAP;
        fn_SMAP = fieldnames(p_SMAP);
        for l = 1:length(fn_SMAP)
            thisFn_SMAP =  fn_SMAP{l};
            switch thisFn_SMAP
                case 'general'
                otherwise
                    h = eval(['h_' thisFn_SMAP]);
                    h.setGuiParameters(p_SMAP.(thisFn_SMAP));
            end
        end

        
    end

    p_LocMoFit = p.(fn{k}).LocMoFit;
    for l = 1:length(p_LocMoFit.parID)
        fitterGUI.fitter.setParArg(p_LocMoFit.parID{k}, 'value', p_LocMoFit.value)
    end
end