function [mask,imff,rangex,rangey]=segmentpointclouds(xy,pixelsize,varargin)
p=parseinput(varargin{:});
x=xy(:,1);y=xy(:,2);
if isempty(p.rangex)
    p.rangex=[min(x)-pixelsize(1) max(x)+pixelsize(1)];
end
if isempty(p.rangey)
    p.rangey=[min(y)-pixelsize(2) max(y)+pixelsize(2)];
end
numpx=(p.rangex(2)-p.rangex(1))/pixelsize(1);
binx=linspace(p.rangex(1),p.rangex(2),numpx);
pixnew(1)=binx(2)-binx(1);

numpy=(p.rangey(2)-p.rangey(1))/pixelsize(2);
biny=linspace(p.rangey(1),p.rangey(2),numpy);
pixnew(2)=biny(2)-biny(1);

im1=histcounts2(x,y,binx,biny);
if p.periodicx
    imh=[im1;im1];
else
    imh=im1;
end
hk=fspecial('gauss',ceil(4*p.sigma),p.sigma);
hkf=fspecial('gauss',ceil(2*p.sigma),p.sigma/2);
if p.sqrt
    imh=sqrt(imh);
end
imf=imfilter(imh,hk,'replicate');
imff=imfilter(imh,hkf,'replicate');

if p.periodicx
    nump=size(im1,1);
    imf=imf(round(nump/4+1:2.5*nump/2),:);
    imff=imff(round(nump/4+1:2.5*nump/2),:);
    im1=imh(round(nump/4+1:2.5*nump/2),:);
end

co1=max(hk(:))*1.1;
imbw=imf>co1;
imbw2=imerode(imbw,strel('disk',ceil(3*p.sigma)));
imrim=imbw;imrim(imbw2)=0;
co2=mean(imff(imrim));
% co2=max(imf(imrim))/2;
imbw3=imff>co2;
mask=imbw3|imbw2;

rangex=p.rangex;rangey=p.rangey;

% imfco=imf;
% imfco(imbw2)=max(hk(:));
% 
% mask=imbw;
% mask=activecontour(imf,imbw2);

% figure(299)
% impl=im1;
% impl(:,:,3)=mask;
% impl(:,:,2)=imrim;
% imagesc(impl);
end

function out=parseinput(varargin)

p = inputParser;
   addParameter(p,'rangex',[],@isnumeric);
   addParameter(p,'rangey',[],@isnumeric);
   addParameter(p,'periodicx',false,@islogical);
   addParameter(p,'sigma',1,@isnumeric);
   addParameter(p,'sqrt',false,@islogical);
   parse(p,varargin{:});
   out=p.Results;
end