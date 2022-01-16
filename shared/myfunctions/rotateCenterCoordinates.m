function [xr,yr]=rotateCenterCoordinates(x,y,time)

c = cov(x-mean(x), y-mean(y));
[a, ev] = eig(c);
[ev,ind] = sort(diag(ev));
[xa, ya] = deal(a(1,ind(end)), a(2,ind(end)));
angle = cart2pol(xa, ya);
[xr,yr]=rotcoord(x-mean(x),y-mean(y),angle);
indx=xr<mean(xr);
if mean(time(indx))>mean(time(~indx)) %increasing position with time
   xr=-xr;
end