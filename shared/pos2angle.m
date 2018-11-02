function angle=pos2angle(pos)

dx=pos(2,1)-pos(1,1);
dy=pos(2,2)-pos(1,2);
angle=atan2(dy,dx)*180/pi;
end