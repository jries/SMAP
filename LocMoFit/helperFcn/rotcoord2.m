function [xo,yo]=rotcoord2(x,y,angle)
    % vectorized rotcoord2
    xo=cos(angle)'.*x+sin(angle)'.*y;
    yo=cos(angle)'.*y-sin(angle)'.*x;
end