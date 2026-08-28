% dc overlay map MINFLUX for scanned beads

locwin=20/2; % This determines how many localizations are used
beaddist=110; % this shoudl be the measured bead distance.

off_est=[10 10];
pixrec=2;

layers=find(g.getPar('sr_layerson'));
locs=g.locData.getloc({'xnm','ynm','time','tid','thi'},'layer',layers,'Position','fov','grouping','ungrouped');

tid=mode(locs.tid(locs.thi==0));
t1=locs.tid==tid & locs.thi==0;
t2=(locs.tid==tid+2 | locs.tid==tid+1) & locs.thi==1;
x1=locs.xnm(t1);y1=locs.ynm(t1);
x2=locs.xnm(t2);y2=locs.ynm(t2);
time1=locs.time(t1);time2=locs.time(t2);
nx=min(min(x2),min(x1)):pixrec:max(max(x2),max(x1));ny=min(min(y2),min(y1)):pixrec:max(max(y2),max(y1));
img1=histcounts2(x1,y1,nx,ny);
imgf1=imgaussfilt(img1,3);
img2=histcounts2(x2,y2,nx,ny);
imgf2=imgaussfilt(img2,3);
clear implot
implot(:,:,1)=imgf1';implot(:,:,2)=imgf2'; implot(:,:,3)=0;
figure(88);hold off
imagesc(implot)
maximaout=maximumfindcall(imgf1');
cutoff=.1;
im=find(maximaout(:,3)>cutoff);
hold on
plot(maximaout(im,1),maximaout(im,2),'w+')


maxx=maximaout(im,1)*pixrec+nx(1);maxy=maximaout(im,2)*pixrec+ny(1);
figure(99);hold off
plot(x1,y1,'r.',x2,y2,'g.')
hold on
plot(maxx,maxy,'k+')
axis ij

searchradius=100;

pos1=zeros(length(im),2);
pos2=zeros(length(im),2);


ims=ceil(sqrt(length(im)));
distmapx=zeros(ims,ims)+NaN;
distmapy=zeros(ims,ims)+NaN;

for k=1:length(im)
    incluster1=(x1-maxx(k)).^2+(y1-maxy(k)).^2 < searchradius;
    incluster2=(x2-maxx(k)).^2+(y2-maxy(k)).^2 < searchradius;
    % figure(89); clf
    % plot(time1(incluster1), y1(incluster1))
    indc1=find(time1>=median(time1(incluster1)),1,'first');
    indc2=find(time2>=median(time1(incluster1)),1,'first');
    pos1(k,1)=mean(x1(indc1-locwin:indc1+locwin));
    pos1(k,2)=mean(y1(indc1-locwin:indc1+locwin));
    pos2(k,1)=mean(x2(indc2-locwin:indc2+locwin));
    pos2(k,2)=mean(y2(indc2-locwin:indc2+locwin));
    bix=round((maxx(k)+beaddist/2)/beaddist);
    biy=round((maxy(k)+beaddist/2)/beaddist);
    bix(bix<1)=1; biy(biy<1)=1;
    distmapx(bix,biy)=pos2(k,1)-pos1(k,1);
    distmapy(bix,biy)=pos2(k,2)-pos1(k,2);
end

nx=isnan(distmapx)|isnan(distmapy);
dxm=mean(distmapx(:),'omitnan');
dym=mean(distmapy(:),'omitnan');
sxm=std(distmapx(:),'omitnan');
sym=std(distmapy(:),'omitnan');
distmapx=distmapx-dxm;distmapy=distmapy-dym; 
% distmapx(nx)=0;distmapy(nx)=0;

figure(91); hold off
imagesc(horzcat(distmapx',distmapy'))
hold on
clear rgbmask
rgbmask(:,:,3)=horzcat(nx',nx');rgbmask(:,:,2)=rgbmask(:,:,3);rgbmask(:,:,1)=rgbmask(:,:,3);
imagesc(rgbmask,"AlphaData",horzcat(nx',nx'))
colorbar
axis equal
title("dx: " + num2str(dxm,'%2.1f') + " nm, dy: " + num2str(dym,'%2.1f') + " nm, " + ...
    "stdx: " + num2str(sxm,'%2.1f') + " nm, stdy: " + num2str(sym,'%2.1f') + " nm")