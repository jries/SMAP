% function countNupsDL
% addSMAPpath
global path 
[f,path]=uigetfile([path filesep '*.*']);
imgr=imageloaderAll([path f]);
% %%

%%
par.cutoffmin=900;
par.mask=0;
par.Rnear=5;
par.sigmaf=0.5;
par.file=f;
kbl=0;
zlen=3;
%%

imgnum=1;
imga=double(imgr.getmanyimages((imgnum-1)*zlen+1:imgnum*zlen,'mat'));
offset=quantile(imga(:),0.02);
img=imga-offset;
%from intensity vs frame fit:
%I(f)=I0*exp(-kbl*f)



for k=1:size(img,3)
    img(:,:,k)=img(:,:,k)/exp(-kbl*(k-1)); %first frame: no bleaching
end

figure
% subplot(2,2,1)
ax1=gca;
% subplot(2,2,2)
% ax2=gca;
% subplot(2,2,3)
% ax0=gca;
ax0=[];ax2=[];
maximaf=getmaxima(img,ax0,ax1,ax2,par);
% colormap(ax0,'gray')

% end
%%
function maximafi=getmaxima(img,ax0,ax1,ax2,par)
cutoffmin=par.cutoffmin;
Rnear=par.Rnear;

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

numbins=100; 
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


%%
maximafi=maxima(indgood,:);

%only take maxima not in first or last frame
% maxp=imghr(maximafi(:,1),maximafi(:,2),1);
[~,indmax]=max(imghr,[],3);
indlin=sub2ind(size(indmax),maximafi(:,1),maximafi(:,2));
indmaxx=indmax(indlin);
mp=size(imghr,3)/2;
% indgood=abs(indmaxx-mp)<=mp-2;
indgood=true(size(indmaxx));
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
hold(ax1,'off')
h=histogram(ax1,mintc,numbins);
% cftool(h.BinEdges(1:end-1)+h.BinWidth/2,h.Values)
xfit=h.BinEdges(1:end-1)+h.BinWidth/2;
yfit=(h.Values);

fitp1=fit(xfit',yfit','gauss1','Robust','LAR');
fitp2=fit(xfit',yfit','gauss2','Robust','LAR','StartPoint',[fitp1.a1,fitp1.b1,fitp1.c1,0,fitp1.b1,fitp1.c1]);
% fitp3=fit(xfit',sqrt(yfit)','gauss3','StartPoint',[fitp.a1,fitp.b1,fitp.c1,fitp.a2*.8,fitp.b2,fitp.c2,fitp.a2*.2,fitp.b1*3,fitp.c2]);
hold(ax1,'on')
plot(ax1,xfit,fitp1(xfit))
plot(ax1,xfit,fitp2(xfit))
legend('hist','single G','double G')

mv=sort([fitp2.b1 fitp2.b2]);
title(ax1,{['maximum at ' num2str(fitp1.b1,4) ', 2G:'  num2str(mv(1),4) ', '  num2str(mv(2),4) ', bg ' num2str(background,2) ', std ' num2str(fitp1.c1/sqrt(2),3)],par.file},...
    'Interpreter','none')
xlabel('intensity at maximum')
ylabel('frequency')
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
