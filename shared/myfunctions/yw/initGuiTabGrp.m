function h = initGuiTabGrp(winName)
    fig = figure('Name',winName);
    h = uitabgroup(fig,'Position',[.01 .01 .98 .98]);
end