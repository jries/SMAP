classdef CME2CSide<interfaces.SEEvaluationProcessor
    properties
        savedevals
    end
    methods
        function obj=CME2CSide(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function makeGui(obj)
            makeGui@interfaces.SEEvaluationProcessor(obj);
%             obj.guihandles.saveimagesb.Callback={@saveimagesb_callback,obj};
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

pard.t1.object=struct('Style','text','String','segment:');
pard.t1.position=[1,1];
pard.t1.Width=1;
pard.segmentl1.object=struct('Style','checkbox','String','layer1','Value',1);
pard.segmentl1.position=[1,2];
pard.segmentl2.object=struct('Style','checkbox','String','layer2','Value',1);
pard.segmentl2.position=[1,3];
pard.segmentl1l2.object=struct('Style','checkbox','String','l1+l2','Value',1);
pard.segmentl1l2.position=[1,4];

pard.t2.object=struct('Style','text','String','layer1');
pard.t2.position=[4,3];
pard.t2.Width=1;
pard.t3.object=struct('Style','text','String','layer2');
pard.t3.position=[4,4];
pard.t3.Width=1;
pard.t4.object=struct('Style','text','String','filter [sx sz]');
pard.t4.position=[5,1];
pard.t4.Width=2;
pard.sigma1.object=struct('Style','edit','String','15 5');
pard.sigma1.position=[5,3];
pard.sigma2.object=struct('Style','edit','String','5 5');
pard.sigma2.position=[5,4];

pard.t5.object=struct('Style','text','String','cutoff factor');
pard.t5.position=[2,1];
pard.t5.Width=2;
pard.cutofffactor.object=struct('Style','edit','String','1');
pard.cutofffactor.position=[2,3.5];

pard.t6.object=struct('Style','text','String','dilation');
pard.t6.position=[3,1];
pard.t6.Width=2;
pard.dilation.object=struct('Style','edit','String','2');
pard.dilation.position=[3,3.5];

pard.plugininfo.type='ROI_Evaluate';

pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi','layer1_','layer2_','se_sitepixelsize'};
end



function out=runintern(obj,p)
if length(p.sigma1)==1
    p.sigma1(2)=p.sigma1(1);
end
if length(p.sigma2)==1
    p.sigma2(2)=p.sigma2(1);
end
roisize=ones(2,1)*p.se_siteroi/2;
roisizeh=p.se_siteroi/2;

[mask,mask1,mask2,mask12,imf1,imf2]=makemaskall(p,obj);
locs1=obj.getLocs({'xnm','ynm','locprecnm','xnmrot','ynmrot'},'layer',1,'size',roisizeh);
locs2=obj.getLocs({'xnm','ynm','locprecnm','xnmrot','ynmrot'},'layer',2,'size',roisizeh);
inmask1=inmask(p,locs1,mask);
inmask2=inmask(p,locs2,mask);

xm1=locs1.xnmrot(inmask1);
ym1=locs1.ynmrot(inmask1);
xm2=locs2.xnmrot(inmask2);
ym2=locs2.ynmrot(inmask2);
%determine parameters

out.N1=sum(inmask1);
out.N2=sum(inmask2);

out.zmed1=median(ym1);out.zmed2=median(ym2);
out.zmedd=out.zmed1-out.zmed2;
out.xmed1=median(xm1);out.xmed2=median(xm2);

out.zm1=mean(ym1);out.zm2=mean(ym2);
out.zmd=out.zmed1-out.zmed2;
out.xm1=mean(xm1);out.xm2=mean(xm2);
sr=p.se_siteroi;
n=-sr/2+p.se_sitepixelsize/2:p.se_sitepixelsize:sr/2;
[X,Y]=meshgrid(n,n);

out.maskmx=sum(sum(mask.*X))/sum(mask(:));
out.maskmz=sum(sum(mask.*Y))/sum(mask(:));

%quantiles
pq=[0.05 0.5 0.95];
if ~isempty(xm1)
    out.qx1=myquantile(xm1,pq);
    out.qz1=myquantile(ym1,pq);
    out.dz1=out.qz1(3)-out.qz1(1);
else
    out.qx1=[0 0 0];
    out.qz1=[0 0 0];
    out.dz1=0;
end
if ~isempty(xm2)
    out.qx2=myquantile(xm2,pq);
    out.qz2=myquantile(ym2,pq);
    out.dz2=out.qz2(3)-out.qz2(1);
else
    out.qx2=[0 0 0];
    out.qz2=[0 0 0];
    out.dz2=0;
end
    out.centdist=out.qz2(2)-out.qz1(2);

%plot all
img1=makeimage(p,locs1.xnmrot,locs1.ynmrot,locs1.locprecnm*p.layer1_.gaussfac);
img2=makeimage(p,locs2.xnmrot,locs2.ynmrot,locs2.locprecnm*p.layer2_.gaussfac);
p.saturation=2;
 hf=obj.setoutput('filtered');
imf=mask2im(p,imf1,imf2,mask1,mask2,mask12,mask);
imagesc([-roisize(1) roisize(1)],[-roisize(2) roisize(2)],imf,'Parent',hf);

 hb=obj.setoutput('image');
imc=mask2im(p,img1,img2,mask1,mask2,mask12,mask);
imagesc([-roisize(1) roisize(1)],[-roisize(2) roisize(2)],imc,'Parent',hb);
%  h2=obj.setoutput('xx');
 h2=hb;
 h2.NextPlot='add';
plot(locs1.xnmrot(inmask1),locs1.ynmrot(inmask1),'bo','Parent',h2);
plot(locs1.xnmrot(~inmask1),locs1.ynmrot(~inmask1),'b.','Parent',h2);
plot(locs2.xnmrot(inmask2),locs2.ynmrot(inmask2),'mo','Parent',h2);
plot(locs2.xnmrot(~inmask2),locs2.ynmrot(~inmask2),'m.','Parent',h2);

plot(out.xmed1,out.zmed1,'ko',out.xmed1,out.zmed1,'wx','Parent',h2);
plot(out.xmed2,out.zmed2,'ko',out.xmed2,out.zmed2,'wx','Parent',h2);
plot(out.xm1,out.zm1,'kd',out.xm1,out.zm1,'w*','Parent',h2);
plot(out.xm2,out.zm2,'kd',out.xm2,out.zm2,'w*','Parent',h2);

plot(out.maskmx,out.maskmz,'wo','Parent',h2);
h2.NextPlot='replace';
axis(h2,'ij');

%make shifted images from filtered localizations
% if ~isnan(out.xmed2)
%     xref=out.xmed2;
% else
%     xref=0;
% end
% if isfield(out,'qz2')
%     yref=out.qz2(3);
% else yref=0;
% end
% 
% x1=locs1.xnmrot(inmask1)-xref;y1=locs1.ynmrot(inmask1)-yref;
% x2=locs2.xnmrot(inmask2)-xref;y2=locs2.ynmrot(inmask2)-yref;
% img1=makeimage(p,x1,y1,locs1.locprecnm(inmask1)*p.layer1_.gaussfac,[],2*roisize);
% img2=makeimage(p,x2,y2,locs2.locprecnm(inmask2)*p.layer2_.gaussfac,[],2*roisize);
% s=size(img1);
% imgt=zeros(s(1),s(2),3);
% imgtp=imgt;
% imgt(:,:,1)=img1;
% imgt(:,:,2)=img2;
% 
% imgtp(:,:,1)=img1/max(img1(:));
% imgtp(:,:,2)=img2/max(img2(:));
% 
% 
%  hc=obj.setoutput('image_shifted');
%  imagesc([-roisize(1) roisize(1)],[-roisize(2) roisize(2)],imgtp,'Parent',hc);
%  out.imageshift=imgt;
end
function imc=mask2im(p,img1,img2,mask1,mask2,mask12,mask)
sim=size(mask1);
imc=zeros(sim(1),sim(2),3);
img1n=img1/max(img1(:))*p.saturation;
img2n=img2/max(img2(:))*p.saturation;
img1n(img1n>1)=1;img2n(img2n>1)=1;

p1=bwperim(mask1);p2=bwperim(mask2);pm=bwperim(mask);p12=bwperim(mask12);
img1n(p1)=1;img2n(p2)=1;
img1n(pm)=1;img2n(pm)=1;
imc(:,:,1)=img1n;
imc(:,:,2)=img2n;
imc(:,:,3)=pm|p12;
imc=imc/max(imc(:));
end
function ind=inmask(p,locs,mask)
roisize=ones(2,1)*p.se_siteroi/2;
pixels=p.se_sitepixelsize;
ind=withinmask(mask',(locs.xnmrot+roisize(1))/pixels,(locs.ynmrot+roisize(2))/pixels);
end
function [mask,mask1,mask2,mask12,img1,img2]=makemaskall(p,obj)
roisize=ones(2,1)*p.se_siteroi/2;
roisizeh=p.se_siteroi/2;
pixels=p.se_sitepixelsize;

locs1=obj.getLocs({'xnm','ynm','locprecnm','xnmrot','ynmrot'},'layer',1,'size',roisizeh,'grouping','ungrouped');
locs2=obj.getLocs({'xnm','ynm','locprecnm','xnmrot','ynmrot'},'layer',2,'size',roisizeh,'grouping','ungrouped');

xm1=locs1.xnmrot;ym1=locs1.ynmrot; % Get coordinates in ROI
xm2=locs2.xnmrot;ym2=locs2.ynmrot; % The same

img1=makeimage(p,xm1,ym1,p.sigma1(1),p.sigma1(2)); % convert coordinates into pixel matrix
img2=makeimage(p,xm2,ym2,p.sigma2(1),p.sigma2(2)); % The same
maxone1=1/(2*pi*p.sigma1(1)*p.sigma1(2))*pixels^2;
maxone2=1/(2*pi*p.sigma2(1)*p.sigma2(2))*pixels^2;

cutoff=p.cutofffactor;
if length(cutoff)==1;
    cutoff(2)=cutoff(1);
end
mask1=makemask(p,img1,maxone1*cutoff(1));



mask=false(size(mask1));
if p.segmentl1
    mask=mask|mask1;
end
if p.segmentl2
    mask2=makemask(p,img2,maxone2*cutoff(2));
    mask=mask|mask2;
else
    mask2=false(size(mask1));
end
if p.segmentl1l2
    mask12=makemask(p,img1/cutoff(1)/maxone1+img2/cutoff(2)/maxone2,1);
    mask=mask|mask12;
else
    mask12=false(size(mask1));
end
end

function [mask,im1,cutoff]=makemask(p,im1,maxone)
% 
% im1=makeimage(p,xm1,ym1,s(1),s(2));

p.take2factor=1.5;
cutoff= maxone;% *p.cutofffactor;
im1bw=im1>cutoff;

% if two largest segments are similar in size
im1bwa1=bwareafilt(im1bw,1);
im1bwa2=bwareafilt(im1bw,2);
if sum(im1bwa2(:))>p.take2factor*sum(im1bwa1(:))
    im1bwa=im1bwa2;
else
    im1bwa=im1bwa1;
end

sel=strel('disk',p.dilation);

im1bwa=imdilate(im1bwa,sel);
mask=imfill(im1bwa,'holes');
end

function im=makeimage(p,xm,ym,sx,sy,roisize)
if nargin<6
    roisize=p.se_siteroi/2;
end
if nargin<5||isempty(sy)
    sy=sx;
end
if length(sx)==1||length(sy)==1
    sx=sx+0*xm; sy=sy+0*xm;
end

pixels=p.se_sitepixelsize;
 range=[-roisize(1) roisize(1)];
 posf.x=xm;posf.y=ym;posf.sx=sx;posf.sy=sy;
im=double(gaussrender_ellipt(posf,range, range, pixels, pixels));
end
