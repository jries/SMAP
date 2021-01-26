h=g.getPar('sr_roihandle');
mask=createMask(h);
pixelsize=g.getPar('sr_pixrec');
area=sum(mask(:)>0)*pixelsize^2;
disp(['area: ' num2str(area/1e6) ' µm^2']);