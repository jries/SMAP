%% Prepare common parameters
% Settings for the simulations
SMAP
p = [];

%% basic settings
p.general.rep = 1;
p.general.saveTo = 'C:\Users\ries\git\SMAP\LocMoFit\docs\images\models\';
p.general.rootPath = 'C:\Users\ries\git\SMAP\LocMoFit\img_model\';
p.general.crop = [0 0]; % from both sides
p.general.view = {'xy'}; % top-view (xy) or side-view (xz, yz)
p.general.scalebar = 50;

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

p.general.thichness = p.settings.se_siteroi; % cross-section

h_guiFile = g.children.guiFile;
%% Get all the listed models
locs.xnm = 0; locs.ynm = 0; locs.znm = 0; locs.locprecnm = 5; locs.locprecznm = 10; locs.layer = 1;
% SE
p_ori = p;
model_imgDefined = img_creation_settings();
for k = 1:length(model_imgDefined)
    oneModel = model_imgDefined{k};
    p_model = img_creation_settings(oneModel);

    fitterGUI = setUpLocMoFit(g, oneModel);
    inp = fitterGUI.getGuiParameters;
    inp.noFit = 1;
    fitterGUI.setGuiParameters(inp);

    p = modPar(fitterGUI, h_SimulateSites, h_settings, p_ori, p_model);
    dimension = fitterGUI.fitter.model{1}.dimension;

    %% Perform simulations
    current_p = h_SimulateSites.getGuiParameters;
    h_SimulateSites.setGuiParameters(current_p);
    h_SimulateSites.setPar('modelType', 'Point');
    for s = 1:p.general.rep
        rng(20220725)
        h_SimulateSites.processgo;
    end

    %% set layer
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

    %% Make plot
    % set ROI
    h = g.getPar('sr_roihandle');
    pos = [p.settings.se_sitefov p.settings.se_sitefov]./1000;
    rot = 90;
    linewidth_roi = p.settings.se_siteroi;
    linelength_roi = p.settings.se_siteroi;
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

            sb = addScalebar(pan(2).axis,'bottom-right', p.settings.se_siteroi./[15 15], 50);
            sb.LineWidth = 3;
        case 3
            numOfView = length(p.general.view);
            fig_pan = figure;
            pan = panel(fig_pan);
            pan.pack('h',numOfView);
            p_general = p.general;
            for j = 1:numOfView
                pan(j).pack('v',2); % one for model, one for simulation
                plotView(pan(j), p.general.view{j}, g, module, fitterGUI, p_general, locs)
            end
            allAx = pan.de.axis;
            hLine = findobj(allAx, 'type', 'line');
            delete(hLine)
            set(allAx, 'YDir', 'normal', 'xtick',[], 'ytick',[],'ytick',[])
            for l = 1:length(allAx)
                allAx(l).YLabel.String = '';
            end
            axis(allAx, 'image')
            axis(allAx, 'on')
            pan.margin = 0;
            pan.margintop = 7;
            pan.marginleft = 7;            
            pan.fontsize = 12;
            pan.de.margin = 0;
            ylabel(pan(1,1).axis,'Model')
            ylabel(pan(1,2).axis,'Simulation')
    end

    %% add the scale bar
    pan_sb = pan(numOfView,2);
    hSb = addScalebar(pan_sb.axis, 'bottom-right', p.settings.se_sitefov*0.05, p.general.scalebar);
    hSb.LineWidth = 3;

    %% Figure size
    fig_pan.Position(3:4) = [7+6+201.5*numOfView 7+201.5*numOfView];

    %% Save
%     dateCreated = ['_c' char(datetime('now','TimeZone','local','Format','yyMMdd'))];
    fName = [p.general.saveTo oneModel '.png'];
    exportgraphics(fig_pan,fName, 'Resolution', 200)

    eval_.removemodule('LocMoFitGUI');
end
%% Functions
function fitterGUI = setUpLocMoFit(g, oneModel)
    se = g.locData.SE;
    eval_ = se.processors.eval;
    eval_.addmodule('LocMoFitGUI');
    idxFitter = find(strcmp('LocMoFitGUI',fieldnames(eval_.children)));
    fitterGUI = eval_.processors{idxFitter};
    fitterGUI.fitter.dataDim = 3;
    fitterGUI.fitter.addModel(functionModel(modelList(oneModel)))
    fitterGUI.fitter.linkedGUI = fitterGUI;
    fitterGUI.updateGUI_fromLocMoFitObj;
end

function p = modPar(fitterGUI, h_SimulateSites, h_settings, p, p_model)
    % 20220818: need to work on this part
    if ~isempty(p_model.SMAP)
        p_model_SMAP = p_model.SMAP;

        % one of 'general','SimulateSites' and 'settings'
        fn_SMAP = fieldnames(p_model_SMAP);
        for l = 1:length(fn_SMAP)
            % update p according to p_model
            thisFn_SMAP =  fn_SMAP{l};
            fn_specificPar = fieldnames(p_model_SMAP.(thisFn_SMAP));
            for k = 1:length(fn_specificPar)
                if isfield(p.(thisFn_SMAP), fn_specificPar{k})
                    p.(thisFn_SMAP).(fn_specificPar{k}) =  p_model_SMAP.(thisFn_SMAP).(fn_specificPar{k});
                else
                    disp(['Invalid parameter name: ' fn_specificPar{k} ' is not defined for ' thisFn_SMAP '.'])
                end
            end
        end   
    end

    % one of 'general','SimulateSites' and 'settings'
    fn = fieldnames(p);
    for l = 1:length(fn)
        thisFn = fn{l};
        switch thisFn
            case 'general'
            otherwise
                h = eval(['h_' thisFn]);
                h.setGuiParameters(p.(thisFn));
        end
    end
    
    %% Set up LocMoFit for simulation
    h_SimulateSites.prepare_LocMoFit(fitterGUI.fitter);
    fitter = h_SimulateSites.getPar('fitter');
    h_SimulateSites.setPar('useDepth', 0);
    h_SimulateSites.setPar('finalROISize',num2str(p.settings.se_siteroi));
    if ~isempty(p_model.LocMoFit)
        p_LocMoFit = p_model.LocMoFit;
        for l = 1:length(p_LocMoFit.parID)
            if ~strcmp(p_LocMoFit.parID{l}, 'm91.sim.numOfMol')
                fitterGUI.fitter.setParArg(p_LocMoFit.parID{l}, 'value', p_LocMoFit.value(l))
            end
            fitter.setParArg(p_LocMoFit.parID{l}, 'value', p_LocMoFit.value(l))
        end
        for l = 1:length(p_LocMoFit.internalSettings)
            fitterGUI.fitter.setModelInternalSetting(1, p_LocMoFit.internalSettings{l}, p_LocMoFit.internalSettings_val(l))
            fitter.setModelInternalSetting(1, p_LocMoFit.internalSettings{l}, p_LocMoFit.internalSettings_val(l))
        end
    end
end

function plotView(pan, view, g, h_viewer3D, fitterGUI, p, locs)
    fitter = fitterGUI.fitter;
    fitter.roiSize = g.getPar('se_siteroi');
    switch view
        case 'xy'
            % model
            ax_mod = fitter.plot(locs, 'pixelSize', 1, 'Projection','xy');
            fig_mod = ax_mod.Parent;
            pan(1).select(ax_mod);
            close(fig_mod);


            % data
            guiPar.theta = 90;
            h_viewer3D.setGuiParameters(guiPar);
            h_viewer3D.run(h_viewer3D.getAllParameters);
            oneView = copy(g.getPar('figurehandle_3Dviewer'));
            pan(2).select(oneView.Children);
            
            close(oneView);
            
            title(pan(1).axis, 'Topview xy')
        case 'xz'
            % model
            ax_mod = fitter.plot(locs, 'pixelSize', 1, 'Projection','xz');
            fig_mod = ax_mod.Parent;
            pan(1).select(ax_mod);
            close(fig_mod);

            % data
            guiPar.theta = 0;
            h_viewer3D.setGuiParameters(guiPar);
            h_viewer3D.run(h_viewer3D.getAllParameters);
            oneView = copy(g.getPar('figurehandle_3Dviewer'));
            pan(2).select(oneView.Children);
            
            close(oneView);
            
            title(pan(1).axis, 'Sideview xz')
        case 'yz'
    end
end
