function dynamicRec_mCME(g, varargin)
    pa = inputParser;
    pa.addParameter('skipRegister',false)
    pa.addParameter('keepCloudsOnly',false)
    pa.parse(varargin{:})
    pa = pa.Results;
    % Disable the layer filtering
    hLayer1 = g.children.guiRender.children.Layer1;
    hLayer1.guihandles.filelistfilter.Value = 0;

    sites = g.locData.SE.sites;
    SMLMModelFit_saveResult = g.getPar('SMLMModelFit_saveResult');
    fit_manager = SMLMModelFit_saveResult.fit_manager;
%     closeAng = fit_manager.getVariable('SMLMModelFitGUI_2.m1.realCloseAngle');
    curvature = fit_manager.getVariable('SMLMModelFitGUI_2.m1.curvature');
    
    lNegCur = curvature<0;
    for k = 1:length(lNegCur)
        if lNegCur(k)
            sites(k).annotation.use = 0;
        end
    end
    
    if pa.skipRegister
        t = 1;
    else
        t = 2;
    end
    
    if pa.keepCloudsOnly
        t = 1;
    end
    
    %% Sorting and registration
    for k = 1:t
        if ~pa.skipRegister
            % enable the first two sorts
            sortROIs = g.children.guiSites.children.Helper.children.SortROIs;
            sortROIs.guihandles.direction1.Value = 2;
            sortROIs.guihandles.sortedit1.String = 'annotation.use';
            sortROIs.guihandles.sortprop1.Value = 3;

            sortROIs.guihandles.direction2.Value = 1;
            sortROIs.guihandles.sortprop2.Value = 4;
            sortROIs.guihandles.sortedit2.String = 'evaluation.SMLMModelFitGUI_2.allParsArg.value(13)';

            % disable all other sorts
            sortROIs.guihandles.sortedit3.String = '';
            sortROIs.guihandles.sortprop3.Value = 1;
            sortROIs.guihandles.sortedit4.String = '';
            sortROIs.guihandles.sortprop4.Value = 1;

            sortROIs.run(sortROIs.getAllParameters);


            SMLMModelFit_saveResult.registerSites_callBack([],[])
        end
        
        %% Dynamic reconstruction
        hDynamicRec = SMLMModelFit_saveResult.dynamicRec_callBack([],[]);
        hDynamicRec.guihandles.spatialTrimXY.String = '70 50';
        hDynamicRec.run;

        hDefineMainCoordinates=plugin('Process','Modify','DefineMainCoordinates');
        p.Vrim=100;

        hDefineMainCoordinates.handle=figure('MenuBar','none','Toolbar','none','Name','DefineMainCoordinates');
        hDefineMainCoordinates.attachPar(g.P);
        hDefineMainCoordinates.attachLocData(g.locData);

        p.Xrim=10;
        hDefineMainCoordinates.setGuiAppearence(p)
        hDefineMainCoordinates.makeGui;
        hDefineMainCoordinates.guihandles.changex.Value = 1;
        hDefineMainCoordinates.guihandles.changey.Value = 1;
        hDefineMainCoordinates.guihandles.changez.Value = 1;
        hDefineMainCoordinates.guihandles.changefield.Value = 1;
        options = hDefineMainCoordinates.guihandles.xfield.String;
        idxField = find(strcmp(options, 'xnmaligned_masterAvgMod'));
        hDefineMainCoordinates.guihandles.xfield.Value = idxField;
        idxField = find(strcmp(options, 'ynmaligned_masterAvgMod'));
        hDefineMainCoordinates.guihandles.yfield.Value = idxField;
        idxField = find(strcmp(options, 'znmaligned_masterAvgMod'));
        hDefineMainCoordinates.guihandles.zfield.Value = idxField;
        hDefineMainCoordinates.run(hDefineMainCoordinates.getAllParameters);

        
        %% Make plots
        
        h = g.getPar('sr_roihandle');
        p = hDynamicRec.getAllParameters;

        hLayer1.guihandles.znm_min.Value = -50;
        hLayer1.guihandles.znm_max.Value = 200;

        linelength_roi = ((p.binNumber+1)*(p.distBetweenBins));
        pos = [-linelength_roi/2+3950 0;linelength_roi/2+3900 0]./1000;

        setPosition(h,pos)

        viewer3DV01 = g.children.guiAnalysis.children.sr3D.processors.Viewer3DV01;

        fig = figure;
        pan = panel(fig);
        pan.pack(2,1)

        viewer3DV01.guihandles.zmax.String = num2str(300);
        viewer3DV01.guihandles.zmin.String = num2str(-100);
        viewer3DV01.guihandles.pixelsize.String = '1,1';

        theta = [90 0];
        for k = 1:length(theta)
            switch theta(k)
                case 90
                    % for
                    g.setPar('linewidth_roi',g.getPar('se_siteroi')-2*p.spatialTrimXY(2));
                case 0
                    g.setPar('linewidth_roi',50);
            end
            viewer3DV01.guihandles.theta.String = num2str(theta(k));

            viewer3DV01.run(viewer3DV01.getAllParameters);
            oneView = copy(g.getPar('figurehandle_3Dviewer'));
            pan(k,1).select(oneView.Children);
            close(oneView);
        end
        ax = pan.de.object;

        set(ax,'XTick',[], 'YTick', [], 'ZTick', []);
        set(ax,'XColor', 'none', 'YColor', 'none', 'ZColor', 'none')
        axis(ax, 'image')
        
        pan.de.margin = 0;
        pan.margin = [2 2 2 2];
        fig.Units = 'centimeters';
        fig.Position(3:4) = 8*[linelength_roi/1000+0.4 1];
    end
end