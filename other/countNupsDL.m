% function countNupsDL
addSMAPpath
global path 
[f,path]=uigetfile([path filesep '*.*']);
imgr=imageloaderAll([path f]);
%%
zlen=5;

imgnum=12
figure(32);
ax0=gca;
figure(33);
for k=1:1
    imgnum=k
imga=double(imgr.getmanyimages((imgnum-1)*zlen+1:imgnum*zlen,'mat'));
offset=0;
img=imga-offset;


ax1=[];
ax2=subplot(5,5,k);
maximaf=getmaxima(img,ax0,ax1,ax2);
end
%

quantile(img(:),0.02)
%%
imgnum=1;
imga=double(imgr.getmanyimages((imgnum-1)*zlen+1:imgnum*zlen,'mat'));
img=imga-offset;
figure(89)
subplot(2,2,1)
ax1=gca;
subplot(2,2,2)
ax2=gca;
subplot(2,2,3)
ax0=gca;
maximaf=getmaxima(img,ax0,ax1,ax2);
colormap(ax0,'gray')

% end
%%
function maximaf=getmaxima(img,ax0,ax1,ax2)
cutoffmin=200;
mask=1;
h=fspecial('gauss',5,0.7);
imgf=imfilter(img,h);
imghr=imresize(imgf,2);
mimgr=max(imghr,[],3);
%%
if 0
figure(39);imagesc(mimgr);
h=imfreehand;
mask=createMask(h);
end
%%
mimg=double(mimgr).*double(mask);
if ~isempty(ax0)
imagesc(ax0,mimg);
end
mmmm=mimg(mimg~=0);
background=quantile(mmmm(:),0.05);
mimg=mimg-background*0;%for now:now bg subtraction
%%
numbins=100; 
maxima=maximumfindcall(mimg);
mint=maxima(:,3);
indgood=mint>cutoffmin;
% cutoffmin=600;
mintc=mint(indgood);

if ~isempty(ax1)
hold(ax1,'off')
h=histogram(ax1,mintc,numbins);
% cftool(h.BinEdges(1:end-1)+h.BinWidth/2,h.Values)
xfit=h.BinEdges(1:end-1)+h.BinWidth/2;
yfit=h.Values;
fitp=fit(xfit',yfit','gauss2','Robust','LAR');
hold(ax1,'on')
plot(ax1,xfit,fitp(xfit))
title(ax1,['maximum at ' num2str(fitp.b1,4) ', bg ' num2str(background) ', std ' num2str(fitp.c1/sqrt(2),3)])
% cftool(h.BinEdges(1:end-1),h.Values)
end
%%
maximaf=maxima(indgood,:);

if ~isempty(ax2)
scatter(ax2,maximaf(:,1),maximaf(:,2),4,maximaf(:,3),'filled')
colormap(ax2,'jet')
axis(ax2,'equal')
axis(ax2,'off')
% ax.CLim=fitp.b1+fitp.c1/sqrt(2)*[-1 1]*2;
ax2.CLim=[500 900];
% colorbar(ax2)
end
end

%%
function statistics
maximac=maxima(mint>cutoffmin,:);
roih=3;
[X,Y]=meshgrid(-roih:roih);
intensity=[];stdx=[];stdy=[];as=[];
for k=size(maximac,1):-1:1
    x=maximac(k,2);y=maximac(k,1);
    if x>roih&& x < size(imghr,2)-roih &&y > roih &&y<size(imghr,1)-roih 
        imsmall=double(mimg(x-roih:x+roih,y-roih:y+roih));
        imsmall=imsmall-min(imsmall(:));
        [as(k),alpha(k)]=asymmetry(imsmall);
        mx=sum(sum(imsmall.*X))/sum(imsmall(:));
        my=sum(sum(imsmall.*Y))/sum(imsmall(:));
        stdx(k)= sum(sum(imsmall.*(X-mx).^2))/sum(imsmall(:));
        stdy(k)= sum(sum(imsmall.*(Y-my).^2))/sum(imsmall(:));
        imrs=imresize(imsmall,8);
        intensity(k)=max(imrs(:));
    end
end



figure(34);
subplot(2,2,1)
plot(as,intensity,'.');
xlabel('asymmetry');
ylabel('intensity');

subplot(2,2,2);
plot(stdx,intensity,'.');

xlabel('stdx');
ylabel('intensity');

subplot(2,2,3);
plot(stdy,intensity,'.');

xlabel('stdy');
ylabel('intensity');

subplot(2,2,4);
plot(stdx./stdy,intensity,'.');

xlabel('stdx/stdy');
ylabel('intensity');

ascutoff=0.15; stdcutoff=3.2;
badind=as>ascutoff|stdx>stdcutoff|stdy>stdcutoff;
badind=badind+mintc'<cutoffmin;

figure(35);
hold off
h=histogram(mintc(~badind),numbins);
% cftool(h.BinEdges(1:end-1),h.Values)

% cftool(h.BinEdges(1:end-1)+h.BinWidth/2,h.Values)
xfit=h.BinEdges(1:end-1)+h.BinWidth/2;
yfit=h.Values;
fitp=fit(xfit',yfit','gauss1','Robust','LAR');
hold on
plot(xfit,fitp(xfit))
title(['Filtered. Maximum at ' num2str(fitp.b1,5) ', bg ' num2str(background)])
end
