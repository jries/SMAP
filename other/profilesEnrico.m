% Calculate radial distribution around beads over time
sites=g.locData.SE.sites;
pixelsize=g.locData.files.file.info.cam_pixelsize_um;
tifffile=g.locData.files.file.tif(1).info.name;
%it got a _1 in the info. strange. Remove that.
tifffile=[tifffile(1:end-4) '.tif'];
tiffreader=mytiffreader(tifffile);
actinstack=double(tiffreader.readall);
actinstack=actinstack-myquantilefast(actinstack,0.02); %BG subtraction

tifffile=g.locData.files.file.tif(2).info.name;
%it got a _1 in the info. strange. Remove that.
tifffile=[tifffile(1:end-4) '.tif'];
tiffreader=mytiffreader(tifffile);
beadstack=double(tiffreader.readall);
beadstack=beadstack-myquantilefast(actinstack,0.02);

%overwrite pixel size:

%% Bead movement
figure(186)
hold off
imagesc(mean(beadstack,3))
hold on
firstframes=3;
lastframes=100;
maxframes=max(g.locData.loc.frame);
for k=length(sites):-1:1
    locs=g.locData.getloc({'xnm','ynm','frame'},'Position',sites(k),'layer',find(g.getPar('sr_layerson')),'grouping','ungrouped');
    ind=find(locs.frame<=firstframes);
    if isempty(ind)
        ind=1;
    end
    xs=mean(locs.xnm(ind));ys=mean(locs.ynm(ind));
    ind=find(locs.frame>=maxframes-lastframes);
    if isempty(ind)
        ind=length(locs.xnm);
    end
    xe=mean(locs.xnm(ind));ye=mean(locs.ynm(ind));
    displ(k,:)=[xe-xs, ye-ys];
    plot([xs,xe]/1000/pixelsize(1),[ys,ye]/1000/pixelsize(2),'d-')
end
dabs=sqrt(displ(:,1).^2+displ(:,2).^2);

figure(185)
histogram(dabs)
%%
%chromatic shift
csx=0;
csy=0;

roisize=30; %2*roisize +1
clear stackcut rsum rsumbin
for k=length(sites):-1:1
    x=round(sites(k).pos(1)/1000/pixelsize(1)+csx);
    y=round(sites(k).pos(2)/1000/pixelsize(2)+csy);

    stackcut(:,:,:,k)=actinstack(y-roisize:y+roisize,x-roisize:x+roisize,:);
%     stackcut(:,:,:,k)=beadstack(y-roisize:y+roisize,x-roisize:x+roisize,:);
    for s=1:size(stackcut,3)
        [rs,norm]=radialsum(stackcut(:,:,s,k));
        rsum(:,s,k)=rs./norm;
    end
end

%%
%mean stack
meanstack=mean(stackcut,4);
fm=figure(187);
imx(meanstack,'Parent',fm)

%%
% use displacement to only average over those

%%
pixelsizeplot=pixelsize(1);
pixelsizeplot=40;

norm =true;
timepoints=10;
df=floor(size(stackcut,3)/timepoints);
figure(188)
subplot(2,2,1);hold off;
col=jet(timepoints+2);col(1:2,:)=[];

radius=((0:size(rsum,1)-1)'+0.5)*pixelsizeplot;
rhd=((0:0.1:size(rsum,1)-1)'+0.5)*pixelsizeplot;
timeaxis=0:timepoints-1;
for t=timepoints:-1:1
    for k=length(sites):-1:1
        rsumbin(:,t,k)=mean(rsum(:,(t-1)*df+1:df*t,k),2);

    end
    rsummean=mean(rsumbin(:,t,:),3);
    if norm %normalize, bleaching
        rsummean=rsummean/rsummean(end);
    end
    stackbin(:,:,t)=mean(meanstack(:,:,(t-1)*df+1:df*t),3);
    subplot(2,2,1)
    plot(radius,rsummean,'Color',col(t,:))
    hold on
    %get HD to work
    rsumhd=interp1(radius,rsummean,rhd,'pchip');
    [imax(t),ind]=max(rsumhd);
    rmax(t)=rhd(ind);
end
title('intensity profile')
xlabel('position (nm)')
ylabel('density (mean intensity)')

subplot(2,2,2)
plot(timeaxis,rmax);
title('radius over time')
xlabel('time')
ylabel('radius (nm)')

subplot(2,2,3)
plot(timeaxis,imax);

title('maximum intensity over time')
xlabel('time')
ylabel('intensity')

fb=figure(187);
imx(stackbin,'Parent',fm)
% stackcutb=beadstack(y-roisize:y+roisize,x-roisize:x+roisize,:);
% function [rs,norm]=callradialsum(img)
% [rs,norm]=radialsum(img);
% end