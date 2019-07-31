classdef CME3DDSpherefit2<interfaces.SEEvaluationProcessor
    properties
        savedevals
        plotfigure
    end
    methods
        function obj=CME3DDSpherefit2(varargin)        
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
pard.text1.object=struct('Style','text','String','slice thickness (nm)');
pard.text1.position=[2,1];
pard.text1.Width=2;
pard.slice.object=struct('Style','edit','String','50');
pard.slice.position=[2,3];
pard.slice.Width=2;

pard.maskon.object=struct('Style','checkbox','String','Use mask. Cutoff:','Value',1);
pard.maskon.position=[3,1];
pard.maskon.Width=2;
pard.maskcutoff.object=struct('Style','edit','String','1');
pard.maskcutoff.position=[3,3];
pard.maskcutoff.Width=2;

pard.textm.object=struct('Style','text','String','Filter mask (pix)');
pard.textm.position=[4,1];
pard.textm.Width=2;
pard.sigmamask.object=struct('Style','edit','String','1');
pard.sigmamask.position=[4,3];
pard.sigmamask.Width=2;

pard.refit.object=struct('Style','checkbox','String','refit sphere','Value',1);
pard.refit.position=[5,1];
pard.refit.Width=2;

pard.saveon.object=struct('Style','checkbox','String','save figure','Value',1);
pard.saveon.position=[6,1];
pard.saveon.Width=2;

pard.zrangecheck.object=struct('Style','checkbox','String','set z center to:','Value',1);
pard.zrangecheck.position=[7,1];
pard.zrangecheck.Width=2;
pard.zrange.object=struct('Style','edit','String','100 ');
pard.zrange.position=[7,3];
pard.zrange.Width=2;

pard.text8.object=struct('Style','text','String','sigma for rendering');
pard.text8.position=[8,1];
pard.text8.Width=2;
pard.renderfilter.object=struct('Style','edit','String','.8');
pard.renderfilter.position=[8,3];
pard.renderfilter.Width=2;


pard.saveimagescheck.object=struct('Style','checkbox','String','');
pard.saveimagescheck.position=[9,1];
pard.saveimagescheck.Width=.3;
pard.saveimagesb.object=struct('Style','pushbutton','String','save as:');
pard.saveimagesb.position=[9,1.3];
pard.saveimagesb.Width=1;
pard.saveimagesdir.object=struct('Style','edit','String','');
pard.saveimagesdir.position=[9,2.4];
pard.saveimagesdir.Width=2.6;
pard.plugininfo.type='ROI_Evaluate';

pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi','layer1_','layer2_','se_sitepixelsize'};
end

function saveimagesb_callback(a,b,obj)
p=obj.getGuiParameters;
f=p.saveimagesdir;
d=uigetdir(f);
if d
obj.guihandles.saveimagesdir.String=d;
end
end

function out=runintern(obj,p)
%use LocIndices for speed if specified
seval= obj.site.evaluation;
if isfield(seval,'LocIndices')
    ingrouped=seval.LocIndices.ingrouped;
    inungrouped=seval.LocIndices.inungrouped;
    addp={'ingrouped',ingrouped,'inungrouped',inungrouped};
else
    addp={};
end

layers=find(p.sr_layerson);
roisize=p.se_siteroi/2;
locs=obj.getLocs({'xnm','ynm','znm'},'layer',1,'size','freeroi',addp{:});

lenbar=20;
ranger=[0 roisize];
rangex=[-roisize roisize];

medianz=median(locs.znm);

rangez=[medianz-roisize medianz+roisize];
if p.zrangecheck
    rangez=[p.zrange-roisize p.zrange+roisize];
end
pixelsx=5;
pixelsz=5;

x=locs.xnm-obj.site.pos(1);
y=locs.ynm-obj.site.pos(2);
z=locs.znm-medianz;

%segment central part
[im,mask]=findInMask(x,y,z,rangex,rangez,p);
if ~p.maskon
    im=true(length(x),1);
    mask=ones(size(mask));
end
 
out.numberOfLocsin3DMask=sum(im);
xn=x(im);
yn=y(im);
zn=z(im);

%determine quantiles.
qx=myquantile(xn,[0.05,0.5,0.95]);
qy=myquantile(yn,[0.05,0.5,0.95]);
qz=myquantile(zn,[0.05,0.5,0.95]);

%positions corrected for median
dx=qx(2);dy=qy(2);
[fi,rn]=cart2pol(xn-dx,yn-dy);
[fi,r]=cart2pol(x-dx,y-dy);

qr=myquantile(rn,0.5);

qrangeall=0.05:0.05:0.95;
out.allqx=myquantile(xn,qrangeall);
out.allqy=myquantile(yn,qrangeall);
out.allqz=myquantile(zn,qrangeall);
out.allqr=myquantile(rn,qrangeall);
out.allrange=qrangeall;
out.qx=qx;
out.qy=qy;
out.qz=qz;
out.qr=qr;




%also for channel 2
if length(layers)>1
locs2=obj.getLocs({'xnm','ynm','znm'},'layer',2,'size','freeroi',addp{:});
x2=locs2.xnm-obj.site.pos(1);
y2=locs2.ynm-obj.site.pos(2);
z2=locs2.znm-medianz;
    mask2=imdilate3(mask,2);
    im2=withinmask(mask2,(x2-rangex(1))/10,(y2-rangex(1))/10,(z2-rangez(1))/10); %pixelsize 10
    xn2=x2(im2);
    yn2=y2(im2);
    zn2=z2(im2);
end

%plot
if obj.display
    immask=makeprojections(xn-dx,yn-dy,zn,rn,p.slice,ranger,rangex,rangez,pixelsx,pixelsz,p.renderfilter,lenbar);
    imout= makeprojections(x,y,z,r,p.slice*4,ranger,rangex,rangez,pixelsx,pixelsz,p.renderfilter,lenbar);
%     h=obj.setoutput('image');
%     hp=h.Parent;
    if isempty(obj.plotfigure)||~isvalid(obj.plotfigure)
        obj.plotfigure=figure;
    end
    hp=obj.plotfigure;
    delete(hp.Children)
    hs1=axes('Parent',hp);subplot(4,6,1,hs1);imagesc(rangex,rangez,immask(6).image,'Parent',hs1); colormap(hs1,hot);
    hs2=axes('Parent',hp);subplot(4,6,2,hs2);imagesc(rangex,rangez,immask(2).image,'Parent',hs2); colormap(hs2,hot);
    hs3=axes('Parent',hp);subplot(4,6,3,hs3);imagesc(rangex,rangez,immask(3).image,'Parent',hs3);colormap(hs3,hot);
    hs4=axes('Parent',hp);subplot(4,6,7,hs4);imagesc(rangex,rangex,immask(1).image,'Parent',hs4); colormap(hs4,hot);
    hs5=axes('Parent',hp);subplot(4,6,8,hs5);imagesc(rangex,rangez,immask(4).image,'Parent',hs5); colormap(hs5,hot);
    hs6=axes('Parent',hp);subplot(4,6,9,hs6);imagesc(rangex,rangez,immask(5).image,'Parent',hs6); colormap(hs6,hot);

    hs4.NextPlot='add';
    rectangle('Position',[qx(1)-dx,qy(1)-dy,qx(3)-qx(1),qy(3)-qy(1)],'EdgeColor',[1 1 1],'Parent',hs4)
    rectangle('Position',[-qr,-qr,2*qr,2*qr],'EdgeColor',[1 1 1],...
        'Curvature',[1,1],...
        'EdgeColor','g','Parent',hs4)
    plot(qx(2)-dx,qy(2)-dy,'c+','Parent',hs4)
    
    hs1.NextPlot='add';
    rectangle('Position',[-qr,qz(1),2*qr,qz(3)-qz(1)],'EdgeColor',[1 1 1],'Parent',hs1)
    plot(0,qz(2),'c+','Parent',hs1)
	
    colormap(hs1,hot)
    sim=size(imout(1).image);
    maskx=bwperim(imresize(squeeze((sum(mask,2))),[sim(1) sim(2)]))';
    masky=bwperim(imresize(squeeze((sum(mask,1))),[sim(1) sim(2)]))';
    maskz=bwperim(imresize(squeeze((sum(mask,3))),[sim(1) sim(2)]))';

    imxy=imout(1).image;imxy(maskz==1)=1;
    imxz=imout(2).image;imxz(maskx==1)=1;
    imyz=imout(3).image;imyz(masky==1)=1;

    ht1=axes('Parent',hp);subplot(4,6,13,ht1);imagesc(rangex,rangez,imout(6).image,'Parent',ht1); colormap(ht1,hot);
    ht2=axes('Parent',hp);subplot(4,6,14,ht2);imagesc(rangex,rangez,imxz,'Parent',ht2);colormap(ht2,hot);
    ht3=axes('Parent',hp);subplot(4,6,15,ht3);imagesc(rangex,rangez,imyz,'Parent',ht3); colormap(ht3,hot);
    ht4=axes('Parent',hp);subplot(4,6,19,ht4);imagesc(rangex,rangex,imxy,'Parent',ht4); colormap(ht4,hot);
    ht5=axes('Parent',hp);subplot(4,6,20,ht5);imagesc(rangex,rangez,imout(4).image,'Parent',ht5); colormap(ht5,hot);
    ht6=axes('Parent',hp);subplot(4,6,21,ht6);imagesc(rangex,rangez,imout(5).image,'Parent',ht6); colormap(ht6,hot);

    hsc=axes('Parent',hp);subplot(3,4,11,hsc);
    hmap2=axes('Parent',hp);subplot(3,4,7,hmap2);
    hmap3=axes('Parent',hp);subplot(3,4,3,hmap3);
    hmap3b=axes('Parent',hp);subplot(3,4,4,hmap3b);
    hmap3c=axes('Parent',hp);subplot(3,4,8,hmap3c);
    hmap3d=axes('Parent',hp);subplot(3,4,12,hmap3d);
    hmap3all={hmap3,hmap3b,hmap3c,hmap3d};
else
    hsc=[];
    hmap2=[];
    hmap3all=[];
end

if p.refit||~isfield(obj.site.evaluation,'CME3DDSpherefit2')
    qr=max(50,abs(qr));
    quantiles=[qx(2),qy(2),qz(2),qr];
    out.spherefit=spherefit(xn,yn,zn,quantiles);
    
    rSphere=out.spherefit(1);
    fitp_sphere=out.spherefit; 
    
%     for plotting the circles in rotated frames:
    xw=fitp_sphere(2)-dx;
    yw=fitp_sphere(3)-dy;
    xppi4=cos(pi/4)*xw+sin(pi/4)*yw;
    xmpi4=cos(-pi/4)*xw+sin(-pi/4)*yw;

    out.plotsphere.xppi4=xppi4;
    out.plotsphere.xmpi4=xmpi4;
    %re-center coordinates
    xc=xn-fitp_sphere(2);
    yc=yn-fitp_sphere(3);
    zc=zn-fitp_sphere(4);

    [tc,pc,rc]=cart2sph(xc,yc,zc);
    out.allqtheta=myquantile(pc,qrangeall);
        
    fitpsc=spherecoverage(pc,hsc);
    out.fitpsc=fitpsc;
    out.fitcoverage_thetatop=min(max(-pi/2,fitpsc.fit(2)),pi/2);
    out.fitcoverage_thetabottom=min(max(-pi/2,fitpsc.fit(3)),pi/2);      
    out.map2D=mapanalysis2D(xc,yc,zc,rSphere,hmap2);
    if length(layers)>1
        xc2=xn2-fitp_sphere(2);
        yc2=yn2-fitp_sphere(3);
        zc2=zn2-fitp_sphere(4);
        out.map3D=mapanalysis3D(xc,yc,zc,rSphere,hmap3all,xc2,yc2,zc2);
    else
        out.map3D=mapanalysis3D(xc,yc,zc,rSphere,hmap3all);
    end

else %use old fit
    oldeval=obj.site.evaluation.CME3DDSpherefit2;
    out=copyfields(oldeval,out);
    out.spherefit=oldeval.spherefit;
    out.allqtheta=oldeval.allqtheta;
    fitp_sphere=oldeval.spherefit;
    fitpsc=oldeval.fitpsc;
    xppi4=oldeval.plotsphere.xppi4;
    xmpi4=oldeval.plotsphere.xmpi4;
    mapanalysis3D(oldeval.map3D,hmap3all);
    spherecoverage(oldeval.fitpsc,hsc);
    mapanalysis2D(oldeval.map2D,hmap2);
end

if obj.display
    hs2.NextPlot='add';
    rectangle('Parent',hs2,'Position',[fitp_sphere(2)-fitp_sphere(1)-dx,fitp_sphere(4)-fitp_sphere(1),2*fitp_sphere(1),2*fitp_sphere(1)],'EdgeColor',[1 1 1], 'Curvature',[1,1])   
    plotwedge(fitp_sphere(2)-dx,fitp_sphere(4),fitp_sphere(1),fitpsc.fit(2),hs2,'c')
    plotwedge(fitp_sphere(2)-dx,fitp_sphere(4),fitp_sphere(1),fitpsc.fit(3),hs2,'c')
    hs3.NextPlot='add';
    rectangle('Parent',hs3,'Position',[fitp_sphere(3)-fitp_sphere(1)-dy,fitp_sphere(4)-fitp_sphere(1),2*fitp_sphere(1),2*fitp_sphere(1)],'EdgeColor',[1 1 1], 'Curvature',[1,1])
    plotwedge(fitp_sphere(2)-dx,fitp_sphere(4),fitp_sphere(1),out.allqtheta(1),hs2,'m')
    plotwedge(fitp_sphere(2)-dx,fitp_sphere(4),fitp_sphere(1),out.allqtheta(end),hs2,'m') 
    hs4.NextPlot='add';
    rectangle('Parent',hs4,'Position',[fitp_sphere(2)-fitp_sphere(1)-dx,fitp_sphere(3)-fitp_sphere(1)-dy,2*fitp_sphere(1),2*fitp_sphere(1)],'EdgeColor',[1 1 1], 'Curvature',[1,1])
    hs5.NextPlot='add';
    rectangle('Parent',hs5,'Position',[xppi4-fitp_sphere(1),fitp_sphere(4)-fitp_sphere(1),2*fitp_sphere(1),2*fitp_sphere(1)],'EdgeColor',[1 1 1], 'Curvature',[1,1])
    hs6.NextPlot='add';
    rectangle('Parent',hs6,'Position',[xmpi4-fitp_sphere(1),fitp_sphere(4)-fitp_sphere(1),2*fitp_sphere(1),2*fitp_sphere(1)],'EdgeColor',[1 1 1], 'Curvature',[1,1])
end

if p.saveimagescheck
    maind=p.saveimagesdir;
    lut=hot(256);
    for k=1:6
        thisd=[maind filesep 'view' num2str(k)];
        if ~exist(thisd,'dir')
            mkdir(thisd)
        end
        fo=[thisd filesep num2str(obj.site.indList) '_S' num2str(obj.site.ID) 'C' num2str(obj.site.info.cell) 'F' num2str(obj.site.info.filenumber) '.tif'];
        imageo=uint8(immask(k).image*255);
        imrgb=ind2rgb(imageo,lut);
        imwrite(imrgb,fo)
    end
end

if p.saveon
    maind=p.saveimagesdir;
        thisd=[maind filesep 'figure' ];
        if ~exist(thisd,'dir')
            mkdir(thisd)
        end
    fo=[thisd filesep num2str(obj.site.indList) '_S' num2str(obj.site.ID) 'C' num2str(obj.site.info.cell) 'F' num2str(obj.site.info.filenumber) '.png'];
    saveas(hp,fo);
end
%filter from layer: or add with filter size and check box
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
if nargin==2
    stat=xc;
    hmap=yc;
else
% nump=256;
% roisize=500; %half size (nm)
% cutofffactor=1.5;
% sigma=2.5;
% 
% rangex=[-roisize roisize];
% rangey=rangex;
% pixelsize=[2*roisize/nump,2*roisize/nump];
pixelsize=[2 2];
sigma=5;
      
% [imbw,imhf]=segmentpointclouds([xc,yc],pixelsize,'rangex',rangex,'rangey',rangey,'periodicx',false,'sigma',sigma);
[imbw,imhf,rangex,rangey]=segmentpointclouds([xc,yc],pixelsize,'sigma',sigma);

% imh=myhist2(xc,yc,pixelsize(1),pixelsize(2),rangex,rangey);
% 
% % nrangex=rangex(1):pixelsize(1):rangex(2);
% % nrangey=rangey(1):pixelsize(2):rangey(2);
% % imhxx=histcounts2(xc,yc,nrangex,nrangey);
% 
% hk=fspecial('gauss',20,sigma);
% imhf=imfilter(sqrt(imh),hk,'replicate');
% 
% cutoffone=max(hk(:))*cutofffactor;
% imbw=imhf>cutoffone;
% imbw=imfill(imbw,'holes');
% 
% %%% try new segmentation
% mask = zeros(size(imhf));
% mask(25:end-25,25:end-25) = 1;
% imbw=activecontour(imhf,mask);


areapix=sum(imbw(:));
areanm=areapix*pixelsize(1)*pixelsize(2);
stat.areanm=areanm;
stat.plot.rangex=rangex;
stat.plot.rangey=rangey;
stat.plot.imhf=imhf;
stat.plot.imbw=imbw;
stat.plot.areanm=areanm;
end
if ~isempty(hmap)
    imagesc(stat.plot.rangex,stat.plot.rangey,(stat.plot.imhf)','Parent',hmap);
    hmap.NextPlot='add';
    colormap(hmap,'jet')
    imagesc(stat.plot.rangex,stat.plot.rangey,(1-stat.plot.imbw')*max(stat.plot.imhf(:)*.5),'AlphaData',0.25,'Parent',hmap);
    hmap.NextPlot='replace';
    title(hmap,['area (nm^2): ' num2str(stat.plot.areanm)]);
end
end

function stat=mapanalysis3D(xc,yc,zc,rSphere,hmap,xc2,yc2,zc2)
if nargin==2 %from previous fit, only plot
    plot3Dmaps(xc,yc)
    return;
end
[tcB,pcB,rcB]=cart2sph(zc,yc,xc);
dualcolor=false;
if nargin>5
    [tcB2,pcB2,rcB2]=cart2sph(zc2,yc2,xc2);
    dualcolor=true;
end
nump=128; %side length of reconstruction in pixels
sigma=3.5;
% cutofffactor=1.1;
filterpixelsize=10;
radiuscutoff=500;
rangex=[-pi pi];
rangey=[-1 1];
pixelsize=[2*pi/nump,2/nump];
scalefactor=1;

medianradius=median(rcB);
sigma=min(3.5,double(max(1,filterpixelsize/medianradius/min(pixelsize))));
isperiodic=true;
refine=false;
[stat,imgs]=mapanalysis3Dint; %(xc,yc,zc,rSphere,hmap,xc2,yc2,zc2,scalefactor);
% drawnow
if stat.rSphere>radiuscutoff
    refine=true;
    isperiodic=false;
    sigma=3;
    pixelsize=[5 5]./stat.rSphere;
%     pos=stat.bwfit.centerpos;
    pos=double([median(tcB) median(sin(pcB))]);
    tcB(tcB>pos(1)+pi)=tcB(tcB>pos(1)+pi)-2*pi;
    tcB(tcB<pos(1)-pi)=tcB(tcB<pos(1)-pi)+2*pi;
    rangex=[pos(1)-nump/2*pixelsize(1) pos(1)+nump/2*pixelsize(1)];
    rangey=[max(-1,pos(2)-nump/2*pixelsize(2)) min(1,pos(2)+nump/2*pixelsize(2))];
    scalefactor=(rangex(2)-rangex(1))*(rangey(2)-rangey(1))/(2*pi*2);
    [stat,imgs]=mapanalysis3Dint;%(xc,yc,zc,rSphere,hmap,xc2,yc2,zc2,scalefactor);
    
end
plot3Dmaps(stat,hmap)

function [stat,imgs]=mapanalysis3Dint%(xc,yc,zc,rSphere,hmap,xc2,yc2,zc2,scalefactor)
%     global imtest
% imh=myhist2(tcB,sin(pcB),pixelsize(1),pixelsize(2),rangex,rangey);
% imh=[imh ;imh];
% hk=fspecial('gauss',ceil(4*sigma),sigma);
% imhf=imfilter(sqrt(imh),hk,'replicate');
% imhf=imhf(nump/4+1:2.5*nump/2,:);
% 
% imtest=imhf;
% 
% cutoffone=max(hk(:))*cutofffactor;
% imbw=imhf>cutoffone;

[imbw,imhf]=segmentpointclouds([tcB,sin(pcB)],pixelsize,'rangex',rangex,'rangey',rangey,'periodicx',isperiodic,'sigma',sigma);

%         coveragepix=sum(imbw(:));
%         coverageFraction=coveragepix/numel(imbw);
statb=areastats(imbw);
statb.coverageFraction=statb.coverageFraction*scalefactor;
stat.coverageFraction=statb.coverageFraction;
stat.coverageArea=stat.coverageFraction*4*pi*rSphere.^2;
stat.rSphere=rSphere;

stat.coordinates.t=tcB;
stat.coordinates.p=pcB;
stat.coordinates.r=rcB;

cbx=statb.centroidx*pixelsize(1)+rangex(1);
cby=statb.centroidy*pixelsize(2)+rangey(1);

statd=areastats(~imbw);
cdx=statd.centroidx*pixelsize(1)+rangex(1);
cdy=statd.centroidy*pixelsize(2)+rangey(1);
if all(imbw==1) %hack: if no hole is found
    cdx=cbx;cdy=cby;
    statd=statb;
end
% statb.main
% statd.main
if statd.main.istop == statb.main.istop
    disp('assignment top bottom not clear')
end

stat.topcoverage=~statd.main.istop;

stat.mainFraction=1-statd.main.areacombined/numel(imbw);
stat.mainArea=stat.mainFraction*4*pi*rSphere.^2;

%determine rotation:smaller area
if statb.main.area<statd.main.area
    mainbw=statb.main.image;
    main=statb.main;
    corrt=0;
else
    mainbw=~statd.main.image;
    main=statd.main;
    corrt=pi;
end
%    corrt=0;
% %     corrp=1;
tfrac=pi-(acos((1-2*stat.mainFraction)));

pp=main.centroidx*pixelsize(1)-pi/2;
pp=main.centroidx*pixelsize(1)+range(1)-pi/2;
% pp=mean(rangex)
tt=asin((main.centroidy*pixelsize(2)-1))+corrt;
startp=double([-tt,pp,tfrac]);
% 
if refine
    if abs(mean(rangex))<1 %somehow to test for negative curvature.
        startp(1)=startp(1)+pi;
    end
%     startp(1)=mean(rangex);
%     startp(2);
end
fitp=fitcoverageangle_bw(mainbw,startp,hmap{2},rangex,rangey);
fitp2=fitcoverageangle_image(sqrt(imhf),fitp,hmap{3},rangex,rangey);

stat.bwfit.theta0=fitp(1);
stat.bwfit.phi0=fitp(2);
stat.bwfit.thetacoverage=fitp(3);
stat.bwfit.fraction=0.5*(1-cos(pi-fitp(3)));
stat.bwfit.centerpos=[cbx cby];
stat.bwfit.centerposd=[cdx cdy];

stat.imfit.theta0=fitp2(1);
stat.imfit.phi0=fitp2(2);
stat.imfit.thetacoverage=fitp2(3);
stat.imfit.fraction=0.5*(1-cos(pi-fitp2(3)));


[ztemp,ycorr]=rotcoord(zc,yc,(fitp2(2)-pi));
[xcorr,zcorr]=rotcoord(xc,ztemp,-(fitp2(1)));


stat.coordinates.xc=xcorr;
stat.coordinates.yc=ycorr;
stat.coordinates.zc=zcorr;
stat.coordinates.x=xc;
stat.coordinates.y=yc;
stat.coordinates.z=zc;
imgs.imhf=imhf;
imgs.imbw=imbw;
imgs.mainbw=mainbw;
if dualcolor
    imh2=myhist2(tcB2,sin(pcB2),pixelsize(1),pixelsize(2),rangex,rangey);
    imh2=[imh2 ;imh2];
    imhf2=imfilter(sqrt(imh2),hk,'replicate');
    imhf2=imhf2(nump/4+1:2.5*nump/2,:);
    [ztemp2,ycorr2]=rotcoord(zc2,yc2,(fitp2(2)-pi));
    [xcorr2,zcorr2]=rotcoord(xc2,ztemp2,-(fitp2(1)));
    stat.coordinates2.xc=xcorr2;
    stat.coordinates2.yc=ycorr2;
    stat.coordinates2.zc=zcorr2;
    stat.coordinates2.x=xc2;
    stat.coordinates2.y=yc2;
    stat.coordinates2.z=zc2;
    imgs.imhf2=imhf2;
    stat.plot.numlocs2=length(xc2);
end
stat.plot.imgs=imgs;
stat.plot.rangex=rangex;
stat.plot.rangey=rangey;
stat.plot.sigma=sigma;
stat.plot.dualcolor=dualcolor;
stat.plot.numlocs=length(xc);
end
end
function plot3Dmaps(stat,hmap)
imgs=stat.plot.imgs;
    imhf=imgs.imhf;
    if ~isempty(hmap) %plot everything
        imagesc(stat.plot.rangex,stat.plot.rangey,(imhf)','Parent',hmap{1});
        hmap{1}.NextPlot='add';
        colormap(hmap{1},'jet')
        imagesc(stat.plot.rangex,stat.plot.rangey,(1-imgs.imbw')*max(imhf(:)*.5),'AlphaData',0.25,'Parent',hmap{1});
        plot(hmap{1},[-pi/2 pi/2],[0 0],'wo')
        plot(hmap{1},stat.bwfit.centerpos(1),stat.bwfit.centerpos(2),'wx')
        plot(hmap{1},stat.bwfit.centerposd(1),stat.bwfit.centerposd(2),'w+')
        hmap{1}.NextPlot='replace';
        title(['coverage area (nm^2): ' num2str(stat.coverageArea,'%6.0f') ', fraction: ' num2str(stat.coverageFraction,2)],'Parent',hmap{1});

        sss=size(imhf);
        outdc=zeros(sss(2),sss(1),3);
        outdc(:,:,1)=imhf'/quantile(imhf(:),.98);
        oneloc=1/stat.plot.sigma^2/2/pi;
        if stat.plot.dualcolor
            imhf2=imgs.imhf2;
            outdc(:,:,2)=imhf2'/max(quantile(imhf2(:),.995),oneloc*1.5);
            numlocs2=stat.locs.numlocs2;
        else
            numlocs2=0;
        end
        outdc(:,:,3)=imgs.mainbw'/2;
        imagesc(stat.plot.rangex,stat.plot.rangey,outdc,'Parent',hmap{4});
        title(hmap{4},['N1: ' num2str(stat.plot.numlocs) ', N2: ' num2str(numlocs2)])   
    end
end


function stat=areastats(imbw)
main=struct('area',0,'areacombined',0,'istop',0,'centroidx',NaN,'centroidy',NaN,'image',[]);
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

% XXSX this does not select the largest one! careful! Above: do for bright
% and dark, and take only smaller of those (truly connected
CC=bwconncomp(imbw);
numpix=cellfun(@numel,CC.PixelIdxList);
[~,idx]=max(numpix);
mainim=false(size(imbw));
mainim(CC.PixelIdxList{idx})=true;
stat.main.image=mainim;

[~,indlargest]=max(stat.area);
dtop=sqrt((stat.centroidx(indlargest)-toppix(1)).^2+(stat.centroidy(indlargest)-toppix(2)).^2);
dbottom=sqrt((stat.centroidx(indlargest)-bottompix(1)).^2+(stat.centroidy(indlargest)-bottompix(2)).^2);
stat.main.istop=dtop<dbottom;
stat.coveragepix=sum(imbw(:));
stat.coverageFraction=stat.coveragepix/numel(imbw);
stat.main.area=stat.area(indlargest);
stat.main.centroidx=stat.centroidx(indlargest);
stat.main.centroidy=stat.centroidy(indlargest);

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
    radius=quantiles(4);
    startptot=[radius,quantiles(1),quantiles(2),quantiles(3)+radius];
    % r x y z
    dz=[0  1000 -1000];
    lb=[0    -200 -200 -2000];
    ub=[2000    200 200 2000];

    tic
    [fitp,r]=implicitfit(fh,startptot,xn,yn,zn,lb,ub,1);
%     indf=1;
    fitbest=fitp;
    for k=1:length(dz)
        startp=startptot; startp(1)=sqrt(startp(1)^2+dz(k)^2);startp(4)=startp(4)+dz(k);
        [fitph,rh]=implicitfit(fh,startp,xn,yn,zn,lb,ub,0);
        if rh<r
%             indf=k;
            fitbest=fitph;
        end
    end
    toc
    tic
%     startp=startptot; startp(1)=sqrt(startp(1)^2+dz(indf)^2); startp(4)=startp(4)+dz(indf);
    [fitp,r]=implicitfit(fh,fitbest,xn,yn,zn,lb,ub,0);
    toc
end

