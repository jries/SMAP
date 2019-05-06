function fitCircleImg(in)
% fit blurred circles to images to extract radii.

if nargin>0 & in==1
    plotresults
    return;
end
global pixelsize results path
results=[];
%fitting exm with circle
 
[file, path]=uigetfile([path  '*.tif']);
file=[path file]
img=imread(file);
img=double(img);
%%

fig=figure(88);hold off
imagesc(img)
dcm_obj = datacursormode(fig);

dcm_obj.UpdateFcn={@fitposition,img,fig};

imf=imfinfo(file);
pixelsize=1000/imf.XResolution %pixelsize in nm. Define yourself if not determined correctly
end

function plotresults
global results pixelsize
rnm=results*pixelsize;
cutoff=3.5;
indb=results(:,1)<cutoff;
figure(100);hold off; histogram(rnm(~indb,1),30);
title([mean(rnm(~indb,1)) std(rnm(~indb,1))]);
xlabel('radius(nm)')
if size(rnm,2)>1

figure(101)
histogram(results(~indb,2),30);
xlabel('Gaussfit sigma (nm)')
end
end
%%
function txt=fitposition(a,cinf,img,fig)
global results
sigmafix=6.8; %SIM
% sigmafix=[];
sigmafix=12; %WF
pos=cinf.Position;
dn=11;
% sroi=2*dn+1;
roi=img(pos(2)-dn:pos(2)+dn,pos(1)-dn:pos(1)+dn);
roi=sqrt(roi-min(roi(:)));
% figure(89);imagesc(roi)
hold on
plot(pos(1),pos(2),'mo')

if 0 % profile
    px=sum(roi,1)';py=sum(roi,2);
    n=(1:sroi)';
    fitp=fit(n,px,'gauss1');
    fitp2=fit(n,py,'gauss1');
    figure(89); plot(n,px,n,fitp(n),n,py,n,fitp2(n));
    title([fitp.c1 fitp2.c1])
    
end
%
startp=[0 0 5 3.3 max(roi(:)) min(roi(:))];
blurrc=(circleblurr(startp,11,sigmafix));
fitc=lsqcurvefit(@circleblurr,startp,dn,roi,[],[],[],sigmafix);
fitroi=circleblurr(fitc,11,sigmafix);
imout=vertcat(horzcat(roi,blurrc),horzcat(fitroi,roi-fitroi));
figure(87); imagesc(imout);
% txt=[num2str(fitc(3),3) ', s=' num2str(fitc(4),3)];
if isempty(sigmafix)
    txt=[num2str(fitc(3),3) ', s=' num2str(fitc(4),3)];
results(end+1,1:2)=fitc(3:4);
else
txt=num2str(fitc(3),3);
results(end+1,1)=fitc(3);
end
end

function out=circleblurr(par,dn,sigmafix)
x=par(1);y=par(2);
R=par(3);sigma=par(4);A=par(5);bg=par(6);
if nargin>2 && ~isempty(sigmafix)
    sigma=sigmafix;
end
n=-dn:dn;
[X,Y]=meshgrid(n,n); 
out=exp(-abs((X-x).^2+(Y-y).^2-(R)^2)/sigma^2)*A+bg;
end