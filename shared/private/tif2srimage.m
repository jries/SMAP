function im=tif2srimage(file,p,form)
if nargin<3
    form='tif';
end
rangext=[p.sr_pos(1)-p.sr_size(1) p.sr_pos(1)+p.sr_size(1)];
rangeyt=[p.sr_pos(2)-p.sr_size(2) p.sr_pos(2)+p.sr_size(2)]; 

rangex=rangext-p.shiftxy_min;
rangey=rangeyt-p.shiftxy_max;
fs=p.ch_filelist.Value;
tnum=p.render_colormode.Value;

%  file=lp.files.file;
 
 fileh=file(fs);
 
 if length(fileh.(form))<tnum||isempty(fileh.(form)(tnum).image)
     im.image=[];%zeros(round(p.sr_sizeRecPix(1)),round(p.sr_sizeRecPix(2)));
     im.rangex=rangex+p.shiftxy_min;
     im.rangey=rangey+p.shiftxy_max;
     return
 end
 
 if strcmp(form,'tif')&&isfield(fileh.(form)(tnum).info,'pixsize') 
    pixsize=fileh.(form)(tnum).info.pixsize;
    roi=fileh.(form)(tnum).info.roi;
 elseif strcmp(form,'tif')&&isfield(fileh.(form)(tnum).info,'cam_pixelsize_um') 
    pixsize=fileh.(form)(tnum).info.cam_pixelsize_um;
    roi=fileh.(form)(tnum).info.roi;
 else
    pixsize=fileh.info.cam_pixelsize_um;
    roi=fileh.info.roi;
%      pixsize(2)=pixsize(1); %XXXXXXX hack to put images together. No idea why. XXXXX
 end

srec=round(p.sr_sizeRecPix);
rangex=rangex+pixsize(1)*1000/2;rangey=rangey+pixsize(end)*1000/2;
rangexpix=rangex/1000/pixsize(1);rangeypix=rangey/1000/pixsize(end);
% rpixrx=[floor(rangexpix(1)) ceil(rangexpix(2))];rpixry=[floor(rangeypix(1)) ceil(rangeypix(2))];
rpixrx=[round(rangexpix(1)) round(rangexpix(2))];rpixry=[round(rangeypix(1)) round(rangeypix(2))];
position=[rpixrx(1),rpixrx(2),rpixry(1),rpixry(2)];

positionc=position;
positionc(1:2)=position(1:2)-roi(1);
positionc(3:4)=position(3:4)-roi(2);

coim=cutoutim(permute(double(fileh.(form)(tnum).image),[2 1 3]),positionc);

magnification=pixsize/p.sr_pixrec*1000;
% magnification=magnification([2 1])
if numel(magnification)>2
    warning('problem in tif2srimage')
end
%%%XXX pixresize: here it is not clear if indices are correct.
% disp('tif2srimage: check rescale line 47')
s=size(coim);
srcoim=imresize(coim,magnification.*s(1:2),'nearest');
% srfinal=imresize(coim,srec,'nearest');

psr=[rangex(:)/p.sr_pixrec-position(1)*magnification(1)+1; rangey(:)/p.sr_pixrec-position(3)*magnification(end)+1];
psr(2)=psr(1)+srec(1)-1;
psr(4)=psr(3)+srec(2)-1;
srfinal=cutoutim(srcoim,psr);

im.image=permute(srfinal,[2,1,3]);
im.rangex=rangext;
im.rangey=rangeyt;