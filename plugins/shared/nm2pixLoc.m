function [loc,locr]=nm2pixLoc(x,y,pixelsize,roi)
loc.xpix=(x/pixelsize(1))-roi(1);
loc.ypix=(y/pixelsize(end))-roi(2);
locr.xpix=round(loc.xpix);
locr.ypix=round(loc.ypix);
end