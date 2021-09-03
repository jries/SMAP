function ax = initGuiTab(h, tabName)
    hTab = findobj(h, 'title', tabName, 'type', 'tab');
    if isempty(hTab)
        hTab = uitab(h,'Title',tabName);
        ax = axes(hTab);
    else
        ax = hTab.Children;
    end
end