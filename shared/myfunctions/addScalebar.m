function h = addScalebar(ax,corner,margin,len)
    xlim = ax.XLim;
    ylim = ax.YLim;
    
    if length(margin) == 1
        margin = repelem(margin,2);
    end
    if isequal(ax.YDir, 'reverse')
        ylim = ylim(end:-1:1);
        margin(2) = -margin(2);
    end
    switch corner
        case 'top-right'
           x = xlim(2)-margin(1);
           x = [x-len x];
           y = ylim(2)-margin(2);
           y = repelem(y,2);
        case 'top-left'
           x = xlim(1)+margin(1);
           x = [x+len x];
           y = ylim(2)-margin(2);
           y = repelem(y,2);
        case 'bottom-right'
           x = xlim(2)-margin(1);
           x = [x-len x];
           y = ylim(1)+margin(2);
           y = repelem(y,2);
        case 'bottom-left'
           x = xlim(1)+margin(1);
           x = [x+len x];
           y = ylim(1)+margin(2);
           y = repelem(y,2);
           
    end
    hold(ax, 'on')
    h = plot(ax, x, y, '-w', 'linewidth',2);
    set(h, 'Tag', 'scale bar')
    hold(ax, 'off')
end