function showROImanager(obj)
processors=obj.locData.SE.processors;
SEpreview=processors.preview;
    if isempty(SEpreview)||~isvalid(SEpreview.handle)
        processors.SEMainGui.make_siteexplorer;
        SEpreview=obj.locData.SE.processors.preview;
    end
        set(SEpreview.handle,'Visible','on')
        figure(SEpreview.handle)