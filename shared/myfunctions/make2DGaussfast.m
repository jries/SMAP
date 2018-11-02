function g=make2DGaussfast(dx,dy,PSFx,PSFy,n)
n=n(:);
lenn=length(n);
gx=exp((-(dx-n).^2)/2/PSFx^2);
gy=exp((-(dy-n).^2)/2/PSFy^2);
g=repmat(gy,[1,lenn]).*repmat(gx',[lenn,1])/pi/2/PSFx/PSFy;
end