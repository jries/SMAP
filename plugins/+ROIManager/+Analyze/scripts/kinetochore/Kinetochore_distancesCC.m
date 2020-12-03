sites=g.locData.SE.sites;
cells=getFieldAsVector(sites,'info.cell');
cellnumbers=unique(cells);

lUsed = getFieldAsVector(sites, 'annotation.use');

cellnumbers_un = [];
%% separate the used and unused sites
if length(sites(~lUsed))>0
    indCells_unused = unique(getFieldAsVector(sites(~lUsed), 'info.cell'));
    cellnumbers_un = [cellnumbers_un cellnumbers(indCells_unused)];
    cellnumbers(indCells_unused) = [];
    %cellnumbers = cellnumbers(cellnumbers~=indCells_unused);
    [imout_un, dx1_un, dy1_un, dx2_un, dy2_un, dx_un, dy_un, indcellAll_un] = extractStat(cellnumbers_un, sites, cells);
end

%% extract statistics
[imout, dx1, dy1, dx2, dy2, dx, dy, indcellAll] = extractStat(cellnumbers, sites, cells);

maxprecision=0.1;
pixrec=sites(1).evaluation.shiftCrossCorrelation.GuiParameters.pixrec;
 imouthr=imresize(imout,pixrec/maxprecision,'cubic');
 [~,linind]=max(imouthr(:));
 [xm,ym]=ind2sub(size(imouthr),linind);
 
dxline=(xm-size(imouthr,1)/2)*maxprecision;
dyline=(ym-size(imouthr,2)/2)*maxprecision;

%% make plots

% Average CC
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

% Per cluster statistics
ID = getFieldAsVector(g.locData.SE.sites(indcellAll),'ID');
if length(sites(~lUsed))>0
    ID_un = getFieldAsVector(g.locData.SE.sites(indcellAll_un),'ID');
end

plotSElink(dx,dy,ID,g.locData.SE,'ok')

if length(sites(~lUsed))>0
    hold on
    plotSElink(dx_un,dy_un,ID_un,g.locData.SE,'xk')
    hold off
end

[stxtx, smx, ssx]=mean2str([dx1;dx2],pr);
[stxty, smy, ssy]=mean2str([dy1;dy2],pr);
title({stxtx; stxty})

subplot(2,2,3)
hold off

if length(sites(~lUsed))>0
    plot([[dx1; dx1_un],[dx2; dx2_un]]',[[dy1; dy1_un],[dy2; dy2_un]]')
else
    plot([dx1 dx2]',[dy1 dy2]')
end
hold on

plot((dx1+dx2)/2,(dy1+dy2)/2,'ko')
if length(sites(~lUsed))>0
    plot((dx1_un+dx2_un)/2,(dy1_un+dy2_un)/2,'kx')
end
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

function [imout, dx1, dy1, dx2, dy2, dx, dy, indcellAll] = extractStat(cellnumbers, sites, cells)
    imout=0*sites(1).evaluation.shiftCrossCorrelation.xcorr;
    dx1=zeros(length(cellnumbers),1);
    dy1=zeros(length(cellnumbers),1);
    dx2=zeros(length(cellnumbers),1);
    dy2=zeros(length(cellnumbers),1);
    for k=length(cellnumbers):-1:1
        indcells=find(cells==cellnumbers(k));
        indcellAll(2*k-1:2*k) = indcells;
        imout=sites(indcells(1)).evaluation.shiftCrossCorrelation.xcorr+imout;
        imout=sites(indcells(2)).evaluation.shiftCrossCorrelation.xcorr+imout;
        % dx and dy for the both clusters in one cell
        dx1(k)=sites(indcells(1)).evaluation.shiftCrossCorrelation.dxline;
        dy1(k)=sites(indcells(1)).evaluation.shiftCrossCorrelation.dyline;
        dx2(k)=sites(indcells(2)).evaluation.shiftCrossCorrelation.dxline;
        dy2(k)=sites(indcells(2)).evaluation.shiftCrossCorrelation.dyline;
    end
    dx=zeros(length(indcellAll),1);
    dy=zeros(length(indcellAll),1);

    for k=1:length(indcellAll)
        dx(k)=sites(indcellAll(k)).evaluation.shiftCrossCorrelation.dxline;
        dy(k)=sites(indcellAll(k)).evaluation.shiftCrossCorrelation.dyline;
    end
end