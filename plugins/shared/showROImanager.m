function showROImanager(obj)
processors=obj.locData.SE.processors;
% SEpreview=processors.preview;
    if ~isfield(processors,'SEpreview') || isempty(SEpreview)||~isvalid(SEpreview.handle)
        processors.SEMainGui.make_siteexplorer;
        
    end
    SEpreview=obj.locData.SE.processors.preview;
        set(SEpreview.handle,'Visible','on')
        fig = figure(SEpreview.handle);
        % Yu-Le edited:
        figPos = fig.Position;
        screesize = get(0,'screensize');
        if screesize(4)<figPos(4)
            fig.Units = "normalized";
            fig.Position(3:4) = [0.2812 0.6250];
        end
        