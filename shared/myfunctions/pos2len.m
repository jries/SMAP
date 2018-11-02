function [len,angle]=pos2len(pos)
dd=pos(2,:)-pos(1,:);
angle=atan2(dd(2),dd(1));
len=sqrt(sum(dd.^2));

