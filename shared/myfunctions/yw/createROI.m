function createROI(g,varargin)
    guiFormat = g.getPar('guiFormat');
    guiFormat.guihandles.roi3.Callback{1}([],[],guiFormat,4,[0 1; 0 0]);
end