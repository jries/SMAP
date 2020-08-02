d=3;
mindistance=11;
minintensity=1e6;
actincutoff=2e5;

%%
if ~exist('path','var')
path='/Volumes/LaCie/Enrico/New Airyscan/mid/';
end
[file, path]=uigetfile([path '*.czi']);


fn=[path file];
%%
imgl=imageloaderOME(fn);
imgl.reader.setSeries(0);
imgsall=imgl.getmanyimages(1:70,'mat');
img1=double(imgsall(:,:,1:2:end));
img2=double(imgsall(:,:,2:2:end));
nim=size(img1,3);

%%
clear ssd
for k=1:nim
    imh=img1(:,:,k);
    ssd(k)=std(imh(:));
end


[~,focalplane]=max(ssd);
figure(88);
subplot(2,2,1)
plot(ssd)
xlabel('z plane')
ylabel('standard deviation')
title('determine focal plane');

img1c=img1(:,:,focalplane-d:focalplane+d);
img2c=img2(:,:,focalplane-d:focalplane+d);
imfit=mean(img1c,3);
img2av=mean(img2c,3);
pfit=struct('roisperfit',100,'iterations',30,'mindistance',mindistance);
tic
disp('start fitting')
locs=fitimagesimple(imfit,pfit);
disp('fitting done');
toc


figure(88);
subplot(2,2,2)
histogram(locs.phot,20);
hold on
plot(minintensity,0,'r*')
hold off
xlabel('photon counts bead')

title('determine cutoff photons');


inphot=locs.phot>minintensity;





%%
%
inx=locs.xpix>roisize+1 & locs.xpix< size(img2av,2)-roisize;
iny=locs.ypix>roisize+1 & locs.ypix< size(img2av,1)-roisize;

use=inphot&inx&iny;


x=round(locs.xpix(use));
y=round(locs.ypix(use));

csx=0;
csy=0;

roisize=30; %2*roisize +1
clear stackcut rsum rsumbin actinsignal
for k=length(x):-1:1

    stackcut(:,:,k)=img2av(y(k)-roisize:y(k)+roisize,x(k)-roisize:x(k)+roisize);
    imh=stackcut(:,:,k);
%     stackcut(:,:,k)=imfit(y(k)-roisize:y(k)+roisize,x(k)-roisize:x(k)+roisize);
        [rs,norm]=radialsum(stackcut(:,:,k));
        rsum(:,k)=rs./norm;
        actinsignal(k)=sum(imh(:));
end



incell=actinsignal>actincutoff;

figure(88)
subplot(2,2,3);
histogram(actinsignal,20)
xlabel('actin signal')
title('determine in/out cell cutoff');
hold on
plot(actincutoff,0,'r*')
hold off

figure(89)
imagesc(imfit);
colorbar
hold on
plot(x,y,'ro');
plot(x(incell),y(incell),'k+');
hold off

figure(90)
imagesc(img2av);
colorbar
hold on
plot(x(incell),y(incell),'ro');
hold off

figure(91);
plot(rsum(:,incell))
hold on
plot(mean(rsum(:,incell),2),'k','LineWidth',3);
hold off

