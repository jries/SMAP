function [xo,yo]=rotcoorddeg(x,y,angle)
xo=cosd(angle)*x+sind(angle)*y;
yo=cosd(angle)*y-sind(angle)*x;