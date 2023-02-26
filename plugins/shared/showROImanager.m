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
        
        if screesize(4)<1000
            factor = [640 800]./screesize(3:4);
            if factor(2)>1
                disp('A larger screen is required for a better user experience.')
            end
            fig.Units = "normalized";
            fig.Position(2:4) = [0.07 factor];
        end
        