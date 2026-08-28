function [gcc,r,g11,g22,info]=paircrosscorrelationfunction(x1,y1,x2,y2,options)
%PAIRCROSSCORRELATIONFUNCTION auto- and cross-pair correlation g(r) of two colors
%
%   [gcc,r,g11,g22,info]=paircrosscorrelationfunction(x1,y1,x2,y2,...)
%
%   Both channels are binned onto the same grid, the correlations are
%   calculated by FFT and normalized by the auto-correlation of the mask
%   (Veatch et al., PLoS ONE 2012). g=1 for a random distribution.
%
%   Input (coordinates in nm):
%       x1,y1   coordinates of channel 1
%       x2,y2   coordinates of channel 2 (optional, then only g11 is returned)
%   Name-value pairs:
%       roihandle   handle to a user defined ROI (imroi) used as mask.
%                   If empty or unusable, the bounding box of the data is used.
%       dr          bin size in nm (default 2)
%       rmax        maximum radius in nm (default: half the field of view)
%       edge        zero padding around the data in nm to avoid wrap-around
%                   of the FFT (default: rmax)
%   Output:
%       gcc     cross-correlation g12(r)
%       r       radius in nm (r(1)=0)
%       g11     auto-correlation of channel 1
%       g22     auto-correlation of channel 2
%       info    struct with the 2D correlation maps, mask, densities, images
%
%   Note: the r=0 bin of the auto-correlations contains the self-pairs of
%   each localization and is usually discarded when plotting.

arguments
    x1
    y1
    x2 = []
    y2 = []
    options.roihandle = []; %handle to user defined ROI for mask
    options.dr = 2
    options.rmax = []
    options.edge = []
end

x1=x1(:); y1=y1(:); x2=x2(:); y2=y2(:);
twocolor=~isempty(x2);
dr=options.dr;

if ~isnumeric(dr) || ~isscalar(dr) || ~(dr>0)
    error('paircrosscorrelationfunction:binwidth', ...
        'the bin width dr must be a positive number, it is %s',mat2str(dr));
end
if numel(x1)<2 || numel(y1)~=numel(x1) || (twocolor && (numel(x2)<2 || numel(y2)~=numel(x2)))
    error('paircrosscorrelationfunction:localizations', ...
        ['need at least 2 localizations with matching x and y in each channel, ' ...
         'channel 1: %d/%d, channel 2: %d/%d. Is a ROI selected that contains localizations?'], ...
        numel(x1),numel(y1),numel(x2),numel(y2));
end

%grid: common for both channels, padded to avoid wrap-around
xall=[x1;x2]; yall=[y1;y2];
rmax=options.rmax;
if isempty(rmax)
    rmax=min(max(xall)-min(xall),max(yall)-min(yall))/2;
end
edge=options.edge;
if isempty(edge)
    edge=rmax;
end
rx=makeedges(min(xall)-edge,max(xall)+edge,dr);
ry=makeedges(min(yall)-edge,max(yall)+edge,dr);

img1=histcounts2(x1,y1,rx,ry);
if twocolor
    img2=histcounts2(x2,y2,rx,ry);
else
    img2=[];
end

%mask: ROI if possible, otherwise bounding box of the localizations
mask=maskfromroi(options.roihandle,rx,ry,dr);
if isempty(mask)
    mask=double((rx(1:end-1)'>=min(xall) & rx(1:end-1)'<=max(xall)) & ...
                (ry(1:end-1) >=min(yall) & ry(1:end-1) <=max(yall)));
end

%maximum radius in pixels, limited by the size of the images
nmax=min(round(rmax/dr),floor(min(size(img1))/2)-1);
r=(0:nmax)'*dr;

%correlations. ifft2 scales all of them by 1/numel, this cancels in g
F1=fft2(img1);
Fm=fft2(mask);
[cm,npix]=radialsum(real(fftshift(ifft2(abs(Fm).^2))),nmax);
Amask=sum(mask(:));
rho1=sum(img1(:))/Amask;

c11=radialsum(real(fftshift(ifft2(abs(F1).^2))),nmax);
g11=c11./cm/rho1^2;

if twocolor
    F2=fft2(img2);
    rho2=sum(img2(:))/Amask;
    c22=radialsum(real(fftshift(ifft2(abs(F2).^2))),nmax);
    c12=radialsum(real(fftshift(ifft2(F1.*conj(F2)))),nmax);
    g22=c22./cm/rho2^2;
    gcc=c12./cm/(rho1*rho2);
else
    rho2=[]; g22=[]; gcc=[];
end

if nargout>4
    info.mask=mask; info.img1=img1; info.img2=img2;
    info.rho1=rho1/dr^2; info.rho2=rho2/dr^2; %densities in 1/nm^2
    info.rx=rx; info.ry=ry; info.dr=dr;
    info.maskcorr=cm; info.npixels=npix;
end
end


function edges=makeedges(mi,ma,dr)
%bin edges, number of bins is odd so that the fftshifted zero-lag pixel is
%exactly in the center of the image
edges=mi:dr:ma+dr;
if mod(numel(edges),2)==1 %even number of bins
    edges(end+1)=edges(end)+dr;
end
end


function [rs,npix]=radialsum(img,nmax)
%sum over rings of width 1 pixel around the center of the fftshifted image,
%rings 0...nmax. Sums (not averages) are returned so that the normalization
%by the mask correlation can be done ring by ring.
s=size(img);
c=floor(s/2)+1;
[ix,iy]=ndgrid((1:s(1))-c(1),(1:s(2))-c(2));
ind=round(sqrt(ix.^2+iy.^2))+1;
inr=ind<=nmax+1;
rs=accumarray(ind(inr),img(inr),[nmax+1 1]);
npix=accumarray(ind(inr),1,[nmax+1 1]);
end


function mask=maskfromroi(roi,rx,ry,px)
%resample the mask of the ROI (drawn on the reconstructed image) onto the
%grid of the histogram images
mask=[];
if isempty(roi)
    return
end
try
    proi=roi.getPosition;
    m=roi.createMask;
    roixim=find(sum(m,1)>0,1,'first');
    roiyim=find(sum(m,2)>0,1,'first');
    wmaskx=find(sum(m,1)>0,1,'last')-roixim;
    wmasky=find(sum(m,2)>0,1,'last')-roiyim;
    switch class(roi)
        case 'imellipse'
            wroi=proi(3:4);
            proimin=proi(1:2);
        otherwise
            wroi=max(proi)-min(proi);
            proimin=min(proi);
    end
    sr_pixrec=(wroi(1)/wmaskx*1000+wroi(2)/wmasky*1000)/2; %nm per pixel of the mask
    maskrs=imresize(double(m'),sr_pixrec/px,'bilinear');
    dxy=floor((-proimin*1000+[rx(1), ry(1)]+[roixim roiyim]*sr_pixrec)/px);
    mask=maskrs(dxy(1)+1:dxy(1)+numel(rx)-1,dxy(2)+1:dxy(2)+numel(ry)-1);
    mask=double(mask>0.5); %round: a >0 threshold would dilate the mask by one
                           %pixel per side and bias g(r) by 2*px/roisize
catch
    mask=[];
end
end
