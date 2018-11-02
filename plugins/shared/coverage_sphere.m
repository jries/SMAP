function im=coverage_sphere(ax,az,thetamin,Elevation,Azimuth)
az=az+pi/2;
ax=ax+pi/2;
[E,A1]=rotateSphericalCoordinates(Elevation,Azimuth,3,-az);
[E2,A2]=rotateSphericalCoordinates(E,A1,2,ax);
% im=2*((T2thetamin-pi/2)-1);
im=E2-(thetamin-pi/2);
co=pi/64;
im=im/co;
im(im<-1)=-1;im(im>1)=1;
end