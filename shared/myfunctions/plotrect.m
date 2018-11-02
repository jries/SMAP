function plotrect(axis,pos,att)
%pos: diagonal corners
if nargin<3
    att=[0 0 0];
end
if iscell(att)
   line([pos(1) pos(3) pos(3) pos(1) pos(1)],[pos(2) pos(2) pos(4) pos(4) pos(2)],'Parent',axis,att{:}) 
else
line([pos(1) pos(3) pos(3) pos(1) pos(1)],[pos(2) pos(2) pos(4) pos(4) pos(2)],'Parent',axis,'Color',att)
end
% plot([pos(1) pos(3) pos(3) pos(1) pos(1)],[pos(2) pos(2) pos(4) pos(4) pos(2)],att,'Parent',axis)