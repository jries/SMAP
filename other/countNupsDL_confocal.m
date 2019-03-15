% function countNupsDL
% addSMAPpath
global path 
[f,path]=uigetfile([path filesep '*.*']);
imgr=imageloaderAll([path f]);
% %%

%%
par.cutoffmin=100; %minimum value. You can start by using a small value and then look at the histogram to determine this value
par.mask=false; %if true, user can draw a mask, only local maxima within this mask are anlayzed
par.Rnear=5; % minimum distance between NPCs
par.sigmaf=0.5; % initial blurring (pixels)
par.file=f; % directory. You can select this with the previous cell
par.regionfilter=false;%if true: only keep clear maxima
par.overwrite=true; %if true: overwrite histogram, otherwise add to existing histogram
kbl=0;
kbl=0.077; %bleaching between frames (exponential decay with kbl [1/frame])
zlen=3; %maximum intensity projection over zlen images
imgnum=1; %first image to be used
%%
imga=double(imgr.getmanyimages((imgnum-1)*zlen+1:imgnum*zlen,'mat'));
offset=quantile(imga(:),0.02);
img=imga-offset;
%from intensity vs frame fit:
%I(f)=I0*exp(-kbl*f)

for k=1:size(img,3)
    img(:,:,k)=img(:,:,k)/exp(-kbl*(k-1)); %first frame: no bleaching
end

figure(8)
if par.overwrite
hold off
else
    hold on
end
ax1=gca;
ax0=[];ax2=[];
maximaf=getmaxima(img,ax0,ax1,ax2,par);

%%
function maximafi=getmaxima(img,ax0,ax1,ax2,par)
cutoffmin=par.cutoffmin;
% Rnear=par.Rnear;

sigmaf=par.sigmaf;
h2=fspecial('disk',3);h1=fspecial('gauss',size(h2,1),sigmaf);
% h=h1-h2/5;%h=h/sum(h(:));
h=h1;

imgf=imfilter(img,h);
imghr=imresize3(imgf,2);
background=quantile(imghr(:),0.05);
imghr=imghr-background;
mimgr=max(imghr,[],3);
%%
if par.mask
figure(39);imagesc(mimgr);
h=imfreehand;
mask=createMask(h);
else
    mask=1;
end
%%
mimg=double(mimgr).*double(mask);
if ~isempty(ax0)
imagesc(ax0,mimg);
end
% mmmm=mimg(mimg~=0);
% 
% mimg=mimg-background;
%%


if 0
maxima=zeros(0,3);
for k=1:size(imghr,3)
    maxima=vertcat(maxima,maximumfindcall(imghr(:,:,k).*mask));
end
indg=brightestinregion(maxima(:,1),maxima(:,2),maxima(:,3),Rnear);
sum(indg(:))/length(indg)
maxima=maxima(indg,:);
% maxima=maximumfindcall(mimg);
else
    
    maxima=maximumfindcall(mimgr.*mask);
end


mint=maxima(:,3);
indgood=mint>cutoffmin;
maximafi=maxima(indgood,:);
binpos=0:50:max(maximafi);
% investigate local neighbourhood, only keep maxima that are clear maxima
if par.regionfilter
roi=3;
cutf=0.85;
r=-roi:roi;
maskh=ones(length(r));
maskh(2:end-1,2:end-1)=0;
indtruemax=false(size(indgood));
for k=1:size(maximafi,1)
    roih=mimgr(maximafi(k,2)+r,maximafi(k,1)+r);
    if all(roih.*maskh<maximafi(k,3)*cutf)
        indtruemax(k)=true;
    end
    
end
figure(90);imagesc(mimgr);hold on;
plot(maximafi(indtruemax,1),maximafi(indtruemax,2),'kx')
hold off
indgood=indtruemax;
else 
    indgood=true(size(maximafi,1),1);
end
%%


%only take maxima not in first or last frame
% maxp=imghr(maximafi(:,1),maximafi(:,2),1);
[~,indmax]=max(imghr,[],3);
indlin=sub2ind(size(indmax),maximafi(:,1),maximafi(:,2));
indmaxx=indmax(indlin);
mp=size(imghr,3)/2;
% indgood=abs(indmaxx-mp)<=mp-2;
% indgood=true(size(indmaxx));
maximaf=maximafi(indgood,:);
if ~isempty(ax2)
scatter(ax2,maximaf(:,1),maximaf(:,2),4,maximaf(:,3),'filled')
colormap(ax2,'jet')
axis(ax2,'equal')
axis(ax2,'off')
% ax.CLim=fitp.b1+fitp.c1/sqrt(2)*[-1 1]*2;
ax2.CLim=[500 2500];
% colorbar(ax2)
end
% zinterp(imghr,maximaf(:,1),maximaf(:,2))

if ~isempty(ax1)

% cutoffmin=600;
% mintc=mint(indgood);
mintc=maximaf(:,3);
% hold(ax1,'off')
h=histogram(ax1,mintc,binpos,'Normalization','probability');

% cftool(h.BinEdges(1:end-1)+h.BinWidth/2,h.Values)
xfit=h.BinEdges(1:end-1)+h.BinWidth/2;
yfit=(h.Values);

fitp1=fit(xfit',yfit','gauss1','Robust','LAR');
fitp2=fit(xfit',yfit','gauss2','Robust','LAR','StartPoint',[fitp1.a1,fitp1.b1,fitp1.c1,0,fitp1.b1*2,fitp1.c1]);
% fitp3=fit(xfit',sqrt(yfit)','gauss3','StartPoint',[fitp.a1,fitp.b1,fitp.c1,fitp.a2*.8,fitp.b2,fitp.c2,fitp.a2*.2,fitp.b1*3,fitp.c2]);
hold(ax1,'on')
plot(ax1,xfit,fitp1(xfit))
plot(ax1,xfit,fitp2(xfit))


mv=sort([fitp2.b1 fitp2.b2]);
title(ax1,{['maximum at ' num2str(fitp1.b1,4) ', 2G:'  num2str(mv(1),4) ', '  num2str(mv(2),4) ', bg ' num2str(background,2) ', std ' num2str(fitp1.c1/sqrt(2),3) ', sem ' num2str(fitp1.c1/sqrt(2)/sqrt(sum(indgood)),3)],par.file},...
    'Interpreter','none')
% title(ax1,{['maximum at ' num2str(fitp1.b1,4) ', bg ' num2str(background,2) ', std ' num2str(fitp1.c1/sqrt(2),3)],par.file},...
%     'Interpreter','none')
xlabel('intensity at maximum')
ylabel('frequency')
legend('hist','single G','double G')
% cftool(h.BinEdges(1:end-1),h.Values)
end


end


function zinterp(img,x,y)
f=(1:size(img,3))';
p=squeeze(img(x(19),y(19),:));

f=f(1:end);
fitp=fit(f,p(f),'poly2');
pmax=-fitp.p2^2/4/fitp.p1+fitp.p3;
fs=f(1):0.2:f(end);
figure(91);plot(f,p(f),fs,fitp(fs),f,0*f+pmax)

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
