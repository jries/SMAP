function ellipse(x,y,rx,ry,varargin)
posf=[x-rx,y-ry,2*rx,2*ry];
rectangle('Position',posf,'Curvature',[1 1],varargin{:})
end