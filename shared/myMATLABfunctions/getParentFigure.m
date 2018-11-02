function fig = getParentFigure(fig)
% if the object is a figure or figure descendent, return the
% figure. Otherwise return [].
while ~isempty(fig) && ~strcmp('figure', get(fig,'type'))
  fig = get(fig,'parent');
end