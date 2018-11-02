function [xo,yo]=rotcoord(x,y,angle)
xo=cos(angle)*x+sin(angle)*y;
yo=cos(angle)*y-sin(angle)*x;