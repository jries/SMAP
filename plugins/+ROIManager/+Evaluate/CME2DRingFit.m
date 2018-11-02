classdef CME2DRingFit<interfaces.SEEvaluationProcessor
    properties
        savedevals
    end
    methods
        function obj=CME2DRingFit(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function makeGui(obj)
            makeGui@interfaces.SEEvaluationProcessor(obj);
            obj.guihandles.saveimagesb.Callback={@saveimagesb_callback,obj};
        end
        function out=run(obj,p)
            out=runintern(obj,p);
         
        end
        
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end

function pard=guidef
pard.layer.object=struct('Style','popupmenu','String','layer1|layer2|layer3|layer4|layer5|layer6');
pard.layer.position=[1,1];
pard.layer.Width=2;
pard.plugininfo.type='ROI_Evaluate';
% pard.text1.object=struct('Style','text','String','slice thickness (nm)');
% pard.text1.position=[2,1];
% pard.text1.Width=2;
% pard.slice.object=struct('Style','edit','String','50');
% pard.slice.position=[2,3];
% pard.slice.Width=2;
% 
% pard.maskon.object=struct('Style','checkbox','String','Use mask. Cutoff:','Value',1);
% pard.maskon.position=[3,1];
% pard.maskon.Width=2;
% pard.maskcutoff.object=struct('Style','edit','String','1');
% pard.maskcutoff.position=[3,3];
% pard.maskcutoff.Width=2;
% 
% pard.textm.object=struct('Style','text','String','Filter mask (pix)');
% pard.textm.position=[4,1];
% pard.textm.Width=2;
% pard.sigmamask.object=struct('Style','edit','String','1');
% pard.sigmamask.position=[4,3];
% pard.sigmamask.Width=2;
% 
% pard.refit.object=struct('Style','checkbox','String','refit sphere','Value',1);
% pard.refit.position=[5,1];
% pard.refit.Width=2;
% 
% pard.saveon.object=struct('Style','checkbox','String','save figure','Value',1);
% pard.saveon.position=[6,1];
% pard.saveon.Width=2;
% 
% pard.zrangecheck.object=struct('Style','checkbox','String','set z center to:','Value',1);
% pard.zrangecheck.position=[7,1];
% pard.zrangecheck.Width=2;
% pard.zrange.object=struct('Style','edit','String','100 ');
% pard.zrange.position=[7,3];
% pard.zrange.Width=2;
% 
% pard.text8.object=struct('Style','text','String','sigma for rendering');
% pard.text8.position=[8,1];
% pard.text8.Width=2;
% pard.renderfilter.object=struct('Style','edit','String','.8');
% pard.renderfilter.position=[8,3];
% pard.renderfilter.Width=2;
% 
% 
% pard.saveimagescheck.object=struct('Style','checkbox','String','');
% pard.saveimagescheck.position=[10,1];
% pard.saveimagescheck.Width=.3;
% pard.saveimagesb.object=struct('Style','pushbutton','String','save as:');
% pard.saveimagesb.position=[10,1.3];
% pard.saveimagesb.Width=1;
% pard.saveimagesdir.object=struct('Style','edit','String','');
% pard.saveimagesdir.position=[10,2.4];
% pard.saveimagesdir.Width=2.6;
end

function saveimagesb_callback(a,b,obj)
p=obj.getGuiParameters.par;
f=p.saveimagesdir;
d=uigetdir(f);
if d
obj.guihandles.saveimagesdir.String=d;
end
end

function out=runintern(obj,p)
roisize=obj.getPar('se_siteroi');%obj.site.sePar.Settings.siteroi/2;
locs=obj.getLocs({'xnm','ynm'},'layer',p.layer.Value,'size',roisize*2);

xm=locs.xnm-obj.site.pos(1);
ym=locs.ynm-obj.site.pos(2);


fh=@circle_implicit;
start_radius=70;
startp=[start_radius,mean(xm),mean(ym)];

% disp('fit')
lb=[0 -2 -2 ]*1000;
ub=[2 2 2 ]*1000;
[fitp,r]=implicitfit(fh,startp,xm,ym,0*xm,lb,ub);


 ho=obj.setoutput('image');
% subplot(1,2,1)

plot(xm,ym,'.','Parent',ho)
ho.NextPlot='add';
pos = [fitp(2)-fitp(1) fitp(3)-fitp(1) 2*fitp(1) 2*fitp(1)];
rectangle('Position',pos,'Curvature',[1 1],'Parent',ho)
ho.NextPlot='replace';
axis(ho,'equal');
% img=site.image;
% sim=size(img);
% 
% pixrec=site.globpar.pixrec/1000;
% range=[- sim(1)/2*pixrec, sim(1)/2*pixrec];
% 
% subplot(1,2,2)
% imagesc(range,range,img);
% hold on
% pos = [fitp(2)-fitp(1) fitp(3)-fitp(1) 2*fitp(1) 2*fitp(1)];
% rectangle('Position',pos,'Curvature',[1 1],'EdgeColor',[0 1 1])
% hold off
% axis equal
% title(fitp(1));

circfactor=1.5;

% for k=1:length(xm)
%     angle(k)=(180/pi)*cart2pol(xm(k),ym(k));
% end

% 
% figure(109)
% histogram(angle,10);
% hold on
% histogram(angle+360,10);
% hold off

out.r1=fitp(1);
out.x0=fitp(2);
out.y0=fitp(3);
out.Ncirc=sum((xm-fitp(2)).^2+(ym-fitp(3)).^2<fitp(1).^2*circfactor^2);

title(['N=' num2str(out.Ncirc) ', r=' num2str(out.r1)],'Parent',ho)


end

function [im,mask]=findInMask(xn,yn,zn,rangex,rangez,p)
%3D reconstruction for mask - change these if masking doesn't work
%correctly
prxy=10;
prz=10;
% cutoff=1;
cutofffactor=p.maskcutoff;
% p.sigmamask=1; %to parameters????
dilationPixels=1; %integer, dilation. Two times along each direction
rangerz=rangez(1):prz:rangez(2);
% rangerx=rangex(1):prxy:rangex(2);
srz=length(rangerz);
srx=round((rangex(2)-rangex(1))/prxy);
outim=zeros(srx,srx,srz);

    
%     posr.x=0;posr.y=0;
%     posr.s=sigmarec*prxy;
%     [~,~,G]=gaussrender(posr,rangex, rangex, prxy, prxy);
    him=fspecial('gaussian',5,p.sigmamask);
      
for k=1:srz-1
    ig=zn>rangerz(k)&zn<rangerz(k+1);
    posr.x=xn(ig);posr.y=yn(ig);
%     posr.s=0*posr.x+sigmarec*prxy;
%     outim(:,:,k)=gaussrender(posr,rangex, rangex, prxy, prxy,[],[],G);
out1=histrender(posr,rangex, rangex, prxy, prxy);
 outim(:,:,k)=filter2(him,out1');
    
end
% cutoff=graythresh(outim)*2;
cutoff=quantile(outim(:),0.995)*cutofffactor;
imbw=0*outim;imbw(outim>cutoff)=1;
cc=bwconncomp(imbw);

lc=zeros(cc.NumObjects,1);
for k=1:cc.NumObjects
    lc(k)=length(cc.PixelIdxList{k});
end
[~,ix]=max(lc);
imbw2=0*outim;
imbw2(cc.PixelIdxList{ix})=1;
% imbw2=bwareaopen(imbw,3);
se = strel('disk',dilationPixels);
imbw3 = imdilate(imbw2>0,se);
imbw3 = imdilate(permute(imbw3,[2 3 1]),se);
imbw3 = imdilate(permute(imbw3,[2 3 1]),se);
imbw3=permute(imbw3,[2 3 1]);
% imbw3(imbw3>0)=1;
%  look3d(double(imbw3));

im=withinmask(imbw3,(xn-rangex(1))/prxy,(yn-rangex(1))/prxy,(zn-rangez(1))/prz);
mask=imbw3;
% sum(im)
end


function xslice=makeprojection(x,y,z,sliced,rangex,rangez,pixelsx,pixelsz,angle,renderfilter,lenbar)

xn=cos(angle)*x+sin(angle)*y;
yn=cos(angle)*y-sin(angle)*x;

confinex=[-sliced sliced];
ig=yn>confinex(1)& yn<confinex(2);
posp.x=xn(ig);
posp.y=z(ig);
% posp.s=posp.x*0+120;

xslice=histrender(posp,rangex, rangez, pixelsx, pixelsz);
if renderfilter>0
    h=fspecial('gaussian',3*ceil(renderfilter),renderfilter);
    xslice=imfilter(xslice,h);
end

% imagesc(rangex,rangez,xslice')!^!
xslice=xslice/max(xslice(:));

xslice(end-1,end-lenbar-1:end-1)=1;
end

function imout= makeprojections(xn,yn,zn,rn,sliced,ranger,rangex,rangez,pixelsx,pixelsz,renderfilter,lenbar)

posp.x=rn;
posp.y=zn;
posp.s=posp.x*0+120;

srim=histrender(posp,ranger, rangez, pixelsx, pixelsz);
[Y,X]=meshgrid(rangez(1)+pixelsz/2:pixelsz:rangez(2),ranger(1)+pixelsx/2:pixelsx:ranger(2));
% subplot(2,3,1)

imn=srim./(X'+pixelsx/4);
% imn=srim';
imn2=[imn(:,end) imn(:,end:-1:2) imn];
if renderfilter>0
    h=fspecial('gaussian',3*ceil(renderfilter),renderfilter);
    imn2=imfilter(imn2,h);
end
imn2=imn2/max(imn2(:));
imn2(end-1,end-lenbar-1:end-1)=1;
imout(6).image=imn2;
imout(6).angle=0;
imout(6).title='radialAverage_rz';


imout(1).image=makeprojection(xn,zn,yn,10,rangex,rangex,pixelsx,pixelsx,0,renderfilter,lenbar);
imout(1).angle=0;
imout(1).title='xy';

si=size(imout(1).image);
imout(1).image(si(1)/2,si(2)/2)=1;
imout(1).image(end-6:end-2,2:6)=1;


% subplot(2,3,2)
imout(2).image=makeprojection(xn,yn,zn,sliced,rangex,rangez,pixelsx,pixelsz,0,renderfilter,lenbar);
imout(2).angle=0;
imout(2).title='xz';
imout(2).image((end-5),2:12)=1;
% subplot(2,3,3)
imout(3).image=makeprojection(xn,yn,zn,sliced,rangex,rangez,pixelsx,pixelsz,pi/2,renderfilter,lenbar);
imout(3).angle=90;
imout(3).title='xz';
imout(3).image(end-12:end-2,5)=1;

% subplot(2,3,5)
imout(4).image=makeprojection(xn,yn,zn,sliced,rangex,rangez,pixelsx,pixelsz,pi/4,renderfilter,lenbar);
imout(4).angle=45;
imout(4).title='xz';
imout(4).image(end-12,2)=1;imout(4).image(end-11,3)=1;imout(4).image(end-10,4)=1;imout(4).image(end-9,5)=1;
% 
% subplot(2,3,6)
imout(5).image=makeprojection(xn,yn,zn,sliced,rangex,rangez,pixelsx,pixelsz,-pi/4,renderfilter,lenbar);
imout(5).angle=135;
imout(5).title='xz';
imout(5).image(end-12,6)=1;imout(5).image(end-11,5)=1;imout(5).image(end-10,4)=1;imout(5).image(end-9,3)=1;


end

%  
function    plotwedge(midx,midy,radius,angle,handle,appearence)
thetm=max(-pi/2,min(pi/2,angle));
plot([midx midx+radius*cos(thetm)],[midy midy+radius*sin(thetm)],appearence,'Parent',handle) 
plot([midx midx-radius*cos(thetm)],[midy midy+radius*sin(thetm)],appearence,'Parent',handle) 
end

function stat=mapanalysis2D(xc,yc,zc,rSphere,hmap)
nump=256;
roisize=500; %half size (nm)
cutofffactor=1.5;
sigma=2.5;

rangex=[-roisize roisize];
rangey=rangex;
pixelsize=[2*roisize/nump,2*roisize/nump];
      
imh=myhist2(xc,yc,pixelsize(1),pixelsize(2),rangex,rangey);

hk=fspecial('gauss',20,sigma);
imhf=imfilter(sqrt(imh),hk,'replicate');
% imhf=imhf(nump/4+1:2.5*nump/2,:);

cutoffone=max(hk(:))*cutofffactor;
imbw=imhf>cutoffone;
imbw=imfill(imbw,'holes');
areapix=sum(imbw(:));
areanm=areapix*pixelsize(1)*pixelsize(2);
stat.areanm=areanm;

if ~isempty(hmap)
    imagesc(rangex,rangey,(imhf)','Parent',hmap);
    hmap.NextPlot='add';
    colormap(hmap,'jet')
    imagesc(rangex,rangey,(1-imbw')*max(imhf(:)*.5),'AlphaData',0.25,'Parent',hmap);
%     plot(hmap,[-pi/2 pi/2],[0 0],'wo')
%     plot(hmap,cbx,cby,'wx')
%      plot(hmap,cdx,cdy,'w+')
    hmap.NextPlot='replace';
%     title(['coverage area (nm^2): ' num2str(stat.coverageArea,'%6.0f') ', fraction: ' num2str(stat.coverageFraction)],'Parent',hmap);
end
end

function stat=mapanalysis3D(xc,yc,zc,rSphere,hmap)
[tc2,pc2,rc2]=cart2sph(zc,yc,xc);
nump=128; %side length of reconstruction in pixels
sigma=3.5;
cutofffactor=1.1;

rangex=[-pi pi];
rangey=[-1 1+1/nump];
pixelsize=[2*pi/nump,2/nump];
      
imh=myhist2(tc2,sin(pc2),pixelsize(1),pixelsize(2),rangex,rangey);
imh=[imh ;imh];
hk=fspecial('gauss',ceil(4*sigma),sigma);
imhf=imfilter(sqrt(imh),hk,'replicate');
imhf=imhf(nump/4+1:2.5*nump/2,:);

cutoffone=max(hk(:))*cutofffactor;
imbw=imhf>cutoffone;

%         coveragepix=sum(imbw(:));
%         coverageFraction=coveragepix/numel(imbw);
statb=areastats(imbw);
stat.coverageFraction=statb.coverageFraction;
stat.coverageArea=stat.coverageFraction*4*pi*rSphere.^2;
stat.rSphere=rSphere;

cbx=statb.centroidx*pixelsize(1)-pi;
cby=statb.centroidy*pixelsize(2)-1;

statd=areastats(~imbw);
cdx=statd.centroidx*pixelsize(1)-pi;
cdy=statd.centroidy*pixelsize(2)-1;
% statb.main
% statd.main
if statd.main.istop == statb.main.istop
    disp('assignment top bottom not clear')
end

stat.topcoverage=~statd.main.istop;

stat.mainFraction=1-statd.main.areacombined/numel(imbw);
stat.mainArea=stat.mainFraction*4*pi*rSphere.^2;
      
if ~isempty(hmap)
    imagesc(rangex,rangey,(imhf)','Parent',hmap);
    hmap.NextPlot='add';
    colormap(hmap,'jet')
    imagesc(rangex,rangey,(1-imbw')*max(imhf(:)*.5),'AlphaData',0.25,'Parent',hmap);
    plot(hmap,[-pi/2 pi/2],[0 0],'wo')
    plot(hmap,cbx,cby,'wx')
     plot(hmap,cdx,cdy,'w+')
    hmap.NextPlot='replace';
    title(['coverage area (nm^2): ' num2str(stat.coverageArea,'%6.0f') ', fraction: ' num2str(stat.coverageFraction)],'Parent',hmap);
end
end

function stat=areastats(imbw)

main=struct('area',0,'areacombined',0,'istop',0);

stat=struct('coveragepix',0,'coverageFraction',0,'main',main,'centroidx',[],'centroidy',[],'area',[],'valid',false);

sim=size(imbw);
nump=sim(1);
bottompix=[nump/4,nump/2];
toppix=[3*nump/4,nump/2];
      
struc=regionprops(imbw,'Centroid','Extrema','Area');
if ~isempty(struc)
for k=1:length(struc)
    stat.centroidx(k)=struc(k).Centroid(2);stat.centroidy(k)=struc(k).Centroid(1);

    ex=struc(k).Extrema;
    stat.touchedge{k}=[min(ex(:,1))<1,max(ex(:,1))>nump,min(ex(:,2))<1,max(ex(:,2))>nump];
    stat.area(k)=struc(k).Area;
    stat.valid=true;
end
[~,indlargest]=max(stat.area);
dtop=sqrt((stat.centroidx(indlargest)-toppix(1)).^2+(stat.centroidy(indlargest)-toppix(2)).^2);
dbottom=sqrt((stat.centroidx(indlargest)-bottompix(1)).^2+(stat.centroidy(indlargest)-bottompix(2)).^2);
stat.main.istop=dtop<dbottom;
stat.coveragepix=sum(imbw(:));
stat.coverageFraction=stat.coveragepix/numel(imbw);
stat.main.area=stat.area(indlargest);

areacombined=stat.area(indlargest);
touchref=stat.touchedge{indlargest};
for k=1:length(stat.touchedge)
    if ~(k==indlargest)
        tochhere=stat.touchedge{k};
        if (touchref(1)&tochhere(1))||(touchref(2)&tochhere(2))||(touchref(3)&tochhere(4))||(touchref(4)&tochhere(3))
            areacombined=areacombined+stat.area(k);
        end  
    end
end

stat.main.areacombined=areacombined;
end
end


function fitp=spherefit(xn,yn,zn,quantiles)
    fh=@sphere_implicit;
    
    startptot=[quantiles(4),quantiles(1),quantiles(2),quantiles(3)+50];
    
    dz=[150  1000 -1000];
%     lb=[0 -2 -2 -2]*1000*0.25;
%     ub=[8 2 2 2]*1000*0.25;
    lb=[0    -200 -200 -2000];
    ub=[2000    200 200 2000];

    tic
    [fitp,r]=implicitfit(fh,startptot,xn,yn,zn,lb,ub,1);
    indf=1;
    for k=1:length(dz)
        startp=startptot; startp(1)=sqrt(startp(1)^2+dz(k)^2);startp(4)=startp(4)+dz(k);
        [fitph,rh]=implicitfit(fh,startp,xn,yn,zn,lb,ub,1);
        if rh<r
            indf=k;
        end
    end
    toc
    tic
    startp=startptot; startp(1)=sqrt(startp(1)^2+dz(indf)^2);startp(4)=startp(4)+dz(indf);
    [fitp,r]=implicitfit(fh,startp,xn,yn,zn,lb,ub,0);
    toc
%     fitp
end

