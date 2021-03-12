function diffViewsNPC(g, pos, rot, linewidth_roi, linelength_roi)
	h = g.getPar('sr_roihandle');
    
    prepos = [0 -linelength_roi/2;0 linelength_roi/2]./1000;
    [prepos(:,1),prepos(:,2)] = rotcoord(prepos(:,1),prepos(:,2),deg2rad(rot));
    pos = prepos+pos;
    
    setPosition(h,pos)
    g.setPar('linewidth_roi',linewidth_roi);
    
    viewer3DV01 = g.children.guiAnalysis.children.sr3D.processors.Viewer3DV01;
    
    fig = figure;
    pan = panel(fig);
    pan.pack(3,1)
    
    viewer3DV01.guihandles.zmax.String = num2str(125);
    viewer3DV01.guihandles.zmin.String = num2str(-125);
    
    theta = [0 -45 -90];
    for k = 1:length(theta)
        viewer3DV01.guihandles.theta.String = num2str(theta(k));
        
        viewer3DV01.run(viewer3DV01.getAllParameters);
        oneView = copy(g.getPar('figurehandle_3Dviewer'));
        pan(k,1).select(oneView.Children);
        close(oneView);
    end
    fig.Position(3:4) = [960 1100];
end