% Calculate radial distribution around beads over time
sites=g.locData.SE.sites;
pixelsize=g.locData.files.file.info.cam_pixelsize_um;
tifffile=g.locData.files.file.tif(1).info.name;
%it got a _1 in the info. strange. Remove that.
tifffile=[tifffile(1:end-6) '.tif'];
tiffreader=mytiffreader(tifffile);
actinstack=double(tiffreader.readall);
actinstack=actinstack-myquantilefast(actinstack,0.02); %BG subtraction

tifffile=g.locData.files.file.tif(2).info.name;
%it got a _1 in the info. strange. Remove that.
tifffile=[tifffile(1:end-6) '.tif'];
tiffreader=mytiffreader(tifffile);
beadstack=double(tiffreader.readall);
beadstack=beadstack-myquantilefast(actinstack,0.02);
%%
%chromatic shift
csx=0;
csy=0;

roisize=30; %2*roisize +1
timepoints=10;
clear stackcut rsum rsumbin
for k=length(sites):-1:1
    x=round(sites(k).pos(1)/1000/pixelsize(1)+csx);
    y=round(sites(k).pos(2)/1000/pixelsize(2)+csy);

    stackcut(:,:,:,k)=actinstack(y-roisize:y+roisize,x-roisize:x+roisize,:);
%     stackcut(:,:,:,k)=beadstack(y-roisize:y+roisize,x-roisize:x+roisize,:);
    for s=1:size(stackcut,3)
        [rs,norm]=callradialsum(stackcut(:,:,s,k));
        rsum(:,s,k)=rs./norm;
    end
end

%%
%mean stack
meanstack=mean(stackcut,4);
fm=figure(187);
imx(meanstack,'Parent',fm)


%%
norm =true;

df=floor(size(stackcut,3)/timepoints);
figure(188)
col=jet(timepoints+2);col(1:2,:)=[];
hold off
for t=timepoints:-1:1
    for k=length(sites):-1:1
        rsumbin(:,t,k)=mean(rsum(:,(t-1)*df+1:df*t,k),2);

    end
    rsummean=mean(rsumbin(:,t,:),3);
    if norm %normalize, bleaching
        rsummean=rsummean/rsummean(end);
    end
    stackbin(:,:,t)=mean(meanstack(:,:,(t-1)*df+1:df*t),3);
    
    plot(rsummean,'Color',col(t,:))
    hold on
end
fb=figure(187);
imx(stackbin,'Parent',fm)

% stackcutb=beadstack(y-roisize:y+roisize,x-roisize:x+roisize,:);
function [rs,norm]=callradialsum(img)
[rs,norm]=radialsum(img);
end