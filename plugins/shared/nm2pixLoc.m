function [loc,locr]=nm2pixLoc(x,y,pixelsize,roi)
loc.x=(x/pixelsize(1))-roi(1);
loc.y=(y/pixelsize(end))-roi(2);
locr.x=round(loc.x);
locr.y=round(loc.y);
end