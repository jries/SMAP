function pd = fAt3Dj_v2(xc, yc, zc,xsize,ysize,zsize,delta_f,coeff)

xc = max(xc,0);
xc = min(xc,xsize-1);

yc = max(yc,0);
yc = min(yc,ysize-1);

zc = max(zc,0);
zc = min(zc,zsize-1);

temp = coeff(xc+1,yc+1,zc+1,:);
pd=sum(delta_f.*(temp(:)));





