function beadlocs=getBeadLocs(x,y,p)
pos.x=x;
pos.y=y;
% dx=p.camPixSizeNm;
dx=p.cam_pixelsize_nm(1);
rangex=dx*round([min(pos.x)-5*dx max(pos.x)+5*dx]/dx);
rangey=dx*round([min(pos.y)-5*dx max(pos.y)+5*dx]/dx);


im=histrender(pos,rangex,rangey,dx,dx);

sigma=2;
h=fspecial('gaussian',15,sigma);
im1f=imfilter(im,h);
maxima=maximumfindcall(im1f);

minbeadsgauss=8;

co=minbeadsgauss/(sigma^2*pi*2);
% co=0.5;
indg=maxima(:,3)>co;
my=maxima(indg,2);
mx=maxima(indg,1);
mynm=my*dx+rangey(1);
mxnm=mx*dx+rangex(1);
    
    

imagesc(rangex,rangey,im1f);
hold on
plot(rangex,rangey,mxnm,mynm,'wo')
hold off

beadlocs.x=mxnm;
beadlocs.y=mynm;
end
