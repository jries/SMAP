function ellipse(x,y,rx,ry,varargin)
posf=[x-r,y-r,2*rx,2*ry];
rectangle('Position',posf,'Curvature',[1 1],varargin{:})
end