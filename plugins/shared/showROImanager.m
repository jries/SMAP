function showROImanager(obj)
processors=obj.locData.SE.processors;
% SEpreview=processors.preview;
    if ~isfield(processors,'SEpreview') || isempty(SEpreview)||~isvalid(SEpreview.handle)
        processors.SEMainGui.make_siteexplorer;
        
    end
    SEpreview=obj.locData.SE.processors.preview;
        set(SEpreview.handle,'Visible','on')
        figure(SEpreview.handle)