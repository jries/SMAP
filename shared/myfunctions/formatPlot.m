
fontname='Arial';
ax=gca;
ax.FontSize=14;
ax.FontName=fontname;
ax.XAxis.FontName=fontname;
ax.XAxis.FontSize=18;
ax.YAxis.FontName=fontname;
ax.YAxis.FontSize=18;

lines=findobj(ax,'type','Line');
for k=1:length(lines)
    lines(k).LineWidth=1;
end