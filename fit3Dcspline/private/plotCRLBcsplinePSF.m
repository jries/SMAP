function plotCRLBcsplinePSF(cspline,ax)
rois=size(cspline.coeff{1},1)-3;
N=5000;bg=50;x=rois/2;y=rois/2;z=(1:size(cspline.coeff{1},3))';
v1=ones(length(z),1);
coords=[v1*x , v1*y , v1* N, v1*bg, z];
crlb=CalSplineCRLB(cspline.coeff{1}, rois, coords);
px=120;
xe=sqrt(crlb(:,1))*px;
ye=sqrt(crlb(:,2))*px;
ze=sqrt(crlb(:,5))*cspline.dz;
zp=(z-cspline.z0)*cspline.dz;
indpl=abs(zp)<700;
plot(ax,zp(indpl),xe(indpl),zp(indpl),ye(indpl),zp(indpl),ze(indpl))
ylim(ax,[0 quantile(ze((indpl)),0.95)*1.1])
xlim(ax,[-700 700]);
legend('x','y','z','location','north')
title(ax,'localization precision for N=5000, bg=50')
xlabel(ax,'z (nm)')
ylabel(ax,'sqrt(CRLB) in nm')
% grid;set(gca, 'YLim', [3 20])
end