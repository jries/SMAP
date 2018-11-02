function circle(x,y,r,varargin)
if r>0
posf=[x-r,y-r,2*r,2*r];
rectangle('Position',posf,'Curvature',[1 1],varargin{:})
end
end