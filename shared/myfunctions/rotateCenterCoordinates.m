function [xr,yr,angleout]=rotateCenterCoordinates(x,y,time)

c = cov(x-mean(x), y-mean(y));
[a, ev] = eig(c);
[ev,ind] = sort(diag(ev));
[xa, ya] = deal(a(1,ind(end)), a(2,ind(end)));
angle = cart2pol(xa, ya);
[xr,yr]=rotcoord(x-mean(x),y-mean(y),angle);
indx=xr<mean(xr);
angleout=angle;
if mean(time(indx))>mean(time(~indx)) %increasing position with time, rotate by pi
   xr=-xr;
   yr=-yr;
   angleout=angle+pi;
end
angleout=mod(angleout,2*pi);