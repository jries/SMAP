function [xr,yr,angleout,centerpoint]=rotateCenterCoordinates(x,y,time,range,angle,centerpoint)

if nargin>3 && length(range)==2%time range passed on
    range=range+min(time);
    indt=time>range(1)&time<range(2);
    xh=x(indt);yh=y(indt);
else
    xh=x;yh=y;
end

if nargin>=5 %angle passed on, only rotate
else

    c = cov(xh-mean(xh), yh-mean(yh));
    [a, ev] = eig(c);
    [ev,ind] = sort(diag(ev));
    [xa, ya] = deal(a(1,ind(end)), a(2,ind(end)));
    angle = cart2pol(xa, ya);
end
if nargin<6
    centerpoint(1)=mean(xh);centerpoint(2)=mean(yh);
end

[xr,yr]=rotcoord(x-centerpoint(1),y-centerpoint(2),angle);
indx=xr<mean(xr);
angleout=angle;
if mean(time(indx))>mean(time(~indx)) %increasing position with time, rotate by pi
   xr=-xr;
   yr=-yr;
   angleout=angle+pi;
end
angleout=mod(angleout,2*pi);