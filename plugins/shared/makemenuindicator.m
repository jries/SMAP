function makemenuindicator(handle,postxt,shift)
if nargin<3
    shift=[0 0];
end

hm=uicontrol('Parent',handle.Parent,'String','=','Style','pushbutton');

units=handle.Units;
handle.Units='pixels';
posp=handle.Position;
handle.Units=units;

w=12; h=12;

hpos=posp(1)+(posp(3)-w)/2;
hside=-0.5;
vpos=posp(2)+(posp(4)-h)/2;
vside=-0.5;
if any(strfind(postxt,'r'))
    hpos=posp(1)+posp(3)-w;
    hside=-1;
end
if any(strfind(postxt,'l'))
    hpos=posp(1);
    hside=1;
end
if any(strfind(postxt,'t'))
    vpos=posp(2)+posp(4)-h;
    vside=-1;
end
if any(strfind(postxt,'b'))
    vpos=posp(2);
    vside=1;
end
if any(strfind(postxt,'i'))
    vpos=vpos;
    hpos=hpos;
end
if any(strfind(postxt,'o'))
    vpos=vpos-vside*w;
    hpos=hpos-hside*w;
end
hm.Position=[hpos+shift(1) vpos+shift(2) w h];
hm.Units='normalized';
% hm.UIContextMenu=handle.UIContextMenu;
hm.TooltipString='Context menu present: right click';
% hm.Enable='inactive';

% hm.HitTest='off';