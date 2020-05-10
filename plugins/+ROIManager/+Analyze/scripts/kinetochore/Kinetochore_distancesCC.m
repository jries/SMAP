sites=g.locData.SE.sites;
imout=0*sites(1).evaluation.shiftCrossCorrelation_4.xcorr;
dx=zeros(length(sites),1);
dy=zeros(length(sites),1);
for k=1:length(sites)
    imout=sites(k).evaluation.shiftCrossCorrelation_4.xcorr+imout;
    dx(k)=sites(k).evaluation.shiftCrossCorrelation_4.dxline;
    dy(k)=sites(k).evaluation.shiftCrossCorrelation_4.dyline;
end
maxprecision=0.1;
pixrec=sites(1).evaluation.shiftCrossCorrelation_4.GuiParameters.pixrec;
 imouthr=imresize(imout,pixrec/maxprecision,'cubic');
 [~,linind]=max(imouthr(:));
 [xm,ym]=ind2sub(size(imouthr),linind);
 
dxline=(xm-size(imouthr,1)/2)*maxprecision;
dyline=(ym-size(imouthr,2)/2)*maxprecision;


dn=(size(imout,1)-1)/2;
nax=(-dn:dn)*pixrec;
figure(88);
subplot(2,2,1);
hold off
imagesc(nax,nax,imout)
hold on
plot(dyline,dxline,'+')
title([dxline, dyline])
subplot(2,2,2);
plot(dx,dy,'o')
title([mean(dx) mean(dy)])
% 
% h=fspecial('gaussian',10,1);
%  imout=filter2(h,imout);
