function rmTick(ax)
    n_Ax = length(ax);
    for k = 1:n_Ax
        oneRmTick(ax(k));
    end
end

function oneRmTick(ax)
    set(ax, 'box', 'off')
    boxLine = findobj(ax, {'tag', 'xbox'},'-or',{'tag','ybox'});
    delete(boxLine);

    xlim = ax.XLim;
    ylim = ax.YLim;

    ybox = yline(ax, ylim(2));
    xbox = xline(ax, xlim(2));

    ybox.Tag = 'ybox';
    xbox.Tag = 'xbox';

    set(ax, 'xlim', xlim, 'ylim', ylim)
end