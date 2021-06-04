function module = openPlugin(g, pPath)
    % pPath should be like
    % {'ROIManager','Analyze','SMLMModelFit_dynamicRec_mCME' }
    module=plugin(pPath{:});
    p.Vrim=100;

    module.handle=figure('MenuBar','none','Toolbar','none','Name',pPath{end});
    module.attachPar(g.P);
    module.attachLocData(g.locData);

    p.Xrim=10;
    module.setGuiAppearence(p)
    module.makeGui;
end