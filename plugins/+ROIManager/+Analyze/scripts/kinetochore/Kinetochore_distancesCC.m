sites=g.locData.SE.sites;
cells=getFieldAsVector(sites,'info.cell');
cellnumbers=unique(cells);


imout=0*sites(1).evaluation.shiftCrossCorrelation.xcorr;
dx1=zeros(length(cellnumbers),1);
dy1=zeros(length(cellnumbers),1);
dx2=zeros(length(cellnumbers),1);
dy2=zeros(length(cellnumbers),1);

lUsed = getFieldAsVector(sites, 'annotation.use');
indCells_unused = unique(getFieldAsVector(sites(~lUsed), 'info.cell'));

cellnumbers = cellnumbers(cellnumbers~=indCells_unused);

for k=length(cellnumbers):-1:1
    indcells=find(cells==cellnumbers(k));
    indcellAll(2*k-1:2*k) = indcells;
    imout=sites(indcells(1)).evaluation.shiftCrossCorrelation.xcorr+imout;
    imout=sites(indcells(2)).evaluation.shiftCrossCorrelation.xcorr+imout;
    dx1(k)=sites(indcells(1)).evaluation.shiftCrossCorrelation.dxline;
    dy1(k)=sites(indcells(1)).evaluation.shiftCrossCorrelation.dyline;
    dx2(k)=sites(indcells(2)).evaluation.shiftCrossCorrelation.dxline;
    dy2(k)=sites(indcells(2)).evaluation.shiftCrossCorrelation.dyline;
end
dx=zeros(length(indcellAll),1);
dy=zeros(length(indcellAll),1);
for k=indcellAll
    dx(k)=sites(k).evaluation.shiftCrossCorrelation.dxline;
    dy(k)=sites(k).evaluation.shiftCrossCorrelation.dyline;
end


maxprecision=0.1;
pixrec=sites(1).evaluation.shiftCrossCorrelation.GuiParameters.pixrec;
 imouthr=imresize(imout,pixrec/maxprecision,'cubic');
 [~,linind]=max(imouthr(:));
 [xm,ym]=ind2sub(size(imouthr),linind);
 
dxline=(xm-size(imouthr,1)/2)*maxprecision;
dyline=(ym-size(imouthr,2)/2)*maxprecision;

pr='%2.1f';

dn=(size(imout,1)-1)/2;
nax=(-dn:dn)*pixrec;
figure(188);
subplot(2,2,1);
hold off
imagesc(nax,nax,imout)
hold on
plot(dyline,dxline,'+')
title([dxline, dyline])
subplot(2,2,2);

ID = getFieldAsVector(g.locData.SE.sites,'ID');
plotSElink(dx,dy,ID,g.locData.SE,'o')
[stxtx, smx, ssx]=mean2str([dx1;dx2],pr);
[stxty, smy, ssy]=mean2str([dy1;dy2],pr);
title({stxtx; stxty})

subplot(2,2,3)
hold off

plot([dx1,dx2]',[dy1,dy2]')
hold on

plot((dx1+dx2)/2,(dy1+dy2)/2,'ko')

plot(0,0,'k*')
[ctxtx, cmx, csx]=mean2str((dx1+dx2)/2,pr);
[ctxty, cmy, csy]=mean2str((dy1+dy2)/2,pr);

title({ctxtx; ctxty})

heading=sprintf('file \t N \t CC av dx \t  dy \t cluster dx \t std \t dy \t std \t cell x \t std \t cell y \t std');
disp(heading)
tab=sprintf('\t');
outstr=[g.getPar('lastSMLFile') tab num2str(length(sites)) tab num2str(dxline) tab ...
    num2str(dyline) tab smx tab ssx tab smy tab ssy...
    tab cmx tab csx tab cmy tab csy];
clipboard('copy',outstr);
% 
% h=fspecial('gaussian',10,1);
%  imout=filter2(h,imout);

function [out,ms,ss]=mean2str(x,pr)
[m,s]=robustMean(x);
out=[num2str(m,pr) '\pm' num2str(s,pr)];
ms=num2str(m);ss=num2str(s);
end