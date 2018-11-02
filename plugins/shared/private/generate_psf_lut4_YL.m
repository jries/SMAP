function imPsfAvg= generate_psf_lut4_YL(fname,pixSzIn, pixSzOut, nMol,frz0,saveName,boxSz,interpZFactor)
%averages and interpolates PSF in XY
% assumes we dont have to interpolate in Z

% nargin= numel(varargin);
MAXSHIFT = 300;%Max PSF shift (nm) when recentring anything larger generates an error
ctrdWindowSzNm = 1500;
boxSzNm =5000;
framePerZPos = 1;
% interpZFactor=1;
ii=1;
% while ii<= nargin
%     if strcmp(varargin{ii},'MaxShift')
%         MAXSHIFT = varargin{ii+1};
%         ii=ii+2;
%     elseif strcmp(varargin{ii},'CentroidWindowSz')
%         ctrdWindowSzNm= varargin{ii+1};
%         ii=ii+2;
%     elseif strcmp(varargin{ii},'BoxSzNm')
%         boxSzNm= varargin{ii+1};
%         ii=ii+2;
%     elseif strcmp(varargin{ii},'FramePerZPos')
%         framePerZPos= varargin{ii+1};
%         ii=ii+2;
%     elseif strcmp(varargin{ii},'InterpZFactor')
%         interpZFactor= varargin{ii+1};
%         ii=ii+2;
%     else
%         ii=ii+1;
%     end
% end

% boxSz = round(boxSzNm/pixSzIn);

% boxSz = 21;
HLFBOX = round(boxSz/2);

%load, preproc - use imreadstack
imstack = double(imreadstack(fname));
%*********************************************************
% %01
% offset = 398.6;
% conversion = 5;
% emGain = 200;

% %02
% offset = 200;
% conversion = 2.8;
% emGain = 1;
% %03
% offset = 200;
% conversion = 2.8;
% emGain = 1;

% %04
% offset = 398.6;
% conversion = 4.553;
% emGain = 1;

% %05
% offset = 200;
% conversion = 2.8;
% emGain = 1;

% %06
% offset = 398.6;
% conversion = 5;
% emGain = 1;

%07
offset = 398.6;
conversion = 5;
emGain = 200;

adu2phot = conversion/emGain;

imstack = (imstack-offset)*adu2phot;

%if multiple frames per Z, average to increase SNR
if framePerZPos>1
    imstack=zFrameAverage(imstack,framePerZPos);
    frz0 = round(frz0/framePerZPos);
end


imz0 = imstack(:,:,frz0);
imLow = imstack(:,:,1);
imHigh = imstack(:,:,end);

%select psf areas  
[posList] = getbeadpos(imz0,imLow,imHigh,nMol,HLFBOX);

%find the centroid
searchLimIn = round(ctrdWindowSzNm/pixSzIn);
posListCtrd = zeros(size(posList));
for ii = 1:nMol
    [posListCtrd(ii,:)]= getPeakPosCentroid(imz0, posList(ii,:), searchLimIn);
end

%interpolate and crop
cropRatio = pixSzIn/pixSzOut;
imPsfInFocus = cropim(imz0,posListCtrd, HLFBOX, cropRatio);


%coalign the centres by cross correlating each PSF to the 1st one
searchLimOut = round(ctrdWindowSzNm/pixSzOut);
posListCorr = crossCorrPsfPos(imPsfInFocus,posListCtrd,cropRatio,searchLimOut);


%crop upscaled in focus im using updated positions 
imPsfInFocus2 = cropim(imz0,posListCorr, HLFBOX, cropRatio);


%find the baseline using the in focus frame - need to subtract
for ii = 1:nMol
    im = imPsfInFocus2(:,:,ii);
    bg(ii) = median(im(:));
end

%find the brightness scale factor - again with in focus frame
% brightness scale factor is mean of all pixels 3sigma above background
% accounts for variations in brightness psf to psf but less noisy than just maximum
imStd = std(double(imz0(:)));
for ii = 1:nMol
    im = imPsfInFocus2(:,:,ii);
    th = bg(ii)+ imStd*3;
    iBright = im(find(im>=th));
    bAvg(ii) = mean(iBright);
end

%crop upscaled stack in xy using updated positions
imPsfStack = cropstack(imstack,posListCorr,HLFBOX,cropRatio);

%<DEBUG
%iw = cast(imPsfStack{2},'uint16');
%imwrite(iw(:,:,2),'testPsf_2.tif');
%for ii=2:nFr
%    imwrite(iw(:,:,ii),'testPsf_2.tif','WriteMode','append');
%end
%keyboard
%DEBUG>

%normalize each PSF
for ii = 1:nMol
    imPsfStackNorm{ii} = (imPsfStack{ii}-bg(ii))/(bAvg(ii)-bg(ii));
end

%average
nFr = size(imstack,3);
imsz=size(imPsfStack{1},1);
imPsfAvg = zeros(imsz,imsz,nFr);
for ii = 1:nFr
    for jj =1:nMol
        imPsfAvg(:,:,ii) = imPsfAvg(:,:,ii) + imPsfStackNorm{jj}(:,:,ii);
    end
end
%normalize
iMin = min(imPsfAvg(:));
iMax = max(imPsfAvg(:));
imPsfAvg= (imPsfAvg-iMin)/(iMax-iMin);

%shift everything so that the in-focus frame CoM is at the centre pixel
psfCtrZ0 = getPsfCtr(imPsfAvg(:,:,frz0),searchLimOut);
imPsfAvg = recentre(imPsfAvg,psfCtrZ0,pixSzOut, MAXSHIFT);

%interpolate along Z if required
if interpZFactor>1
    imPsfAvgSave = interpZ(imPsfAvg, interpZFactor);
    nFrSave = size(imPsfAvgSave,3);
else
    imPsfAvgSave = imPsfAvg;
    nFrSave = nFr;
end

iw = im2uint16(imPsfAvgSave);
imwrite(iw(:,:,1),saveName);
for ii=2:nFrSave
    imwrite(iw(:,:,ii),saveName,'WriteMode','append');
end



%%
%plot the wobble curve
psfCtr = getPsfCtr(imPsfAvg,searchLimOut);
szY = size(imPsfAvg,1);
ctrY = floor(szY/2)+1;
imCtr = repmat([ctrY,ctrY],nFr,1);
psfShift = psfCtr-imCtr;
psfShiftNm = psfShift*pixSzOut;
d=sqrt(sum(psfShift.^2,2));
d=d*pixSzOut;
figure;plot(d);
xlabel('z frame');ylabel('PSF CoM dist from centre pix nm');
figure;
hold all;
plot(psfShiftNm(:,1),psfShiftNm(:,2));
plot(psfShiftNm(frz0,1),psfShiftNm(frz0,2),'rx');
xlabel('z frame');ylabel('PSF shift from centre pix nm');

%%
% psf = iw(:,:,2:36);
% psf = iw(:,:,5:41);
% psf = iw(:,:,10:35);
%**************************************
% %01
% psf = iw(:,:,33-18:33+18);
% %02
% psf = iw(:,:,22-18:22+18);
% %03
% psf = iw(:,:,23-18:23+18);
% %04
% psf = iw(:,:,22-18:22+18);
% %05
% psf = iw(:,:,19-18:19+18);
% %06
% psf = iw(:,:,15-14:15+14);
% %06_01
% psf = iw(:,:,16-14:16+14);

%07
psf = iw(:,:,33-12:33+13);



pixel_size = 138;
type = '3D';
zvals = -900:50:900;
zmin = -900;
zmax=900;
I = find(saveName=='\',1,'last');
currentFolder = saveName(1:I);
save([currentFolder 'avgpsf_xy' num2str(pixSzOut)],'psf','pixel_size','type','zmin','zmax','zvals');
%-----------------------------------
function [posList]= getbeadpos(imz0,imLow,imHigh,nMol,HLFBOX)

imsz = size(imz0);

h1=figure; 
subplot(1,3,1);
imagesc(imLow);
axis equal;
axis([0 imsz(2) 0 imsz(1)]);
hold all;
subplot(1,3,2);
imagesc(imz0);
axis equal;
axis([0 imsz(2) 0 imsz(1)]);
hold all;
subplot(1,3,3);
imagesc(imHigh);
axis equal;
axis([0 imsz(2) 0 imsz(1)]);
hold all;

h2=figure;
hold all;
imagesc(imz0);
axis equal;
axis([0 imsz(2) 0 imsz(1)]);
set(gca,'YDir','reverse');

for ii = 1:nMol
    figure(h2)
    [x,y] = ginput(1);
    posList(ii,:)=[x,y];
    box{ii} = [x-HLFBOX,y-HLFBOX;...
            x+HLFBOX,y-HLFBOX;...
            x+HLFBOX,y+HLFBOX;...
            x-HLFBOX,y+HLFBOX;...
            x-HLFBOX,y-HLFBOX];

    plot(x,y,'kx');
    plot(box{ii}(:,1),box{ii}(:,2),'k-');
    fprintf('%d of %d\n',ii,nMol); 
end

for ii = 1:nMol
    for jj=1:3
        figure(h1)
        subplot(1,3,jj);
        plot(x,y,'kx');
        plot(box{ii}(:,1),box{ii}(:,2),'k-');
    end
end

%--------------------------------------------------------
function [pos amplitude]= getPeakPosCentroid(im, posGuess, windowRadius)
% crop small subimage around each im
% get centroid of each sub image


[windowLim] = getSubWindow(im,posGuess,windowRadius);
xLim  = windowLim(1:2);
yLim  = windowLim(3:4);

subIm = im(yLim(1):yLim(2),xLim(1):xLim(2));

posGuessSubIm = [posGuess(1) - xLim(1), posGuess(2)-yLim(1)];
[posSubIm, amplitude] = getCentroid(subIm);

pos = posSubIm + [xLim(1),yLim(1)]-[1,1];

%hold off;
%imagesc(subIm);
%hold all;
%plot(posGuessSubIm(:,1),posGuessSubIm(:,2),'kx');
%plot(posSubIm(:,1),posSubIm(:,2),'mx');
%pause
%
%hold off;
%imagesc(im);
%hold all;
%plot(posGuess(:,1),posGuess(:,2),'kx');
%plot(pos(:,1),pos(:,2),'gx');
%keyboard

%--------------------------------------------------
function [windowLim posGuessSubIm] = getSubWindow(im,posGuess,windowRadius);

[sizey sizex] = size(im);
X0=posGuess(1);
Y0=posGuess(2);

%round X0, Y0 to use as matrix locations
X0_int = round(X0); 
Y0_int = round(Y0);
windowRadius = round(windowRadius); %radius should already be an integer anyway

% setup the limits of the cropped image
xstart =  X0_int-windowRadius;
xfinish = X0_int+windowRadius;
ystart =  Y0_int-windowRadius;
yfinish = Y0_int+windowRadius;
% check if any of the limits are out of bounds - if so, skip that point
if (xstart<1)
  xstart=1;
end
%if (xstart > sizex)
%if (xfinish<1) 
if (xfinish > sizex)
  xfinish=sizex;
end
if (ystart<1)
  ystart=1; 
end
%if (ystart > sizey)
%if (yfinish<1)
if (yfinish > sizey) 
  yfinish = sizey;
end
windowLim = [xstart, xfinish,ystart,yfinish];
%---------------------------------------------
function [pos, amplitude] = getCentroid(im);
im = im-min(im(:));%offset the image to avoid baseline "zero" values contributing to the avg

[sizey sizex] = size(im);
totalIntensity = sum(im(:));
posx = sum(sum(im,1).*[1:sizex]) ./ totalIntensity;
posy = sum(sum(im,2).*[1:sizey]') ./ totalIntensity;

pos = [posx,posy];
amplitude = im(round(posy),round(posx));

%----------------------------------------------------
function imPsf = cropim(im,posList, HLFBOX, cropRatio)

%interpolate the image
imBig = imresize(im,cropRatio);
posListBig = posList*cropRatio;
hlfBoxBig = round(HLFBOX*cropRatio);
boxSz = hlfBoxBig*2+1;

nPsf = size(posList,1);
imPsf = zeros(boxSz, boxSz,nPsf);
for ii=1:nPsf
    x = round(posListBig(ii,1));
    y = round(posListBig(ii,2));
    imPsf(:,:,ii) = imBig(y-hlfBoxBig:y+hlfBoxBig, x-hlfBoxBig:x+hlfBoxBig);
end
    
%----------------------------------------------------------------------------------------------
function imPsfStack = cropstack(imstack,posList,HLFBOX,cropRatio);

nFr = size(imstack,3);
nPsf = size(posList,1);

%for each frame crop around the particle centre
for ii = 1:nFr
    iStackTmp{ii} = cropim(imstack(:,:,ii),posList, HLFBOX, cropRatio);
end

%rearrange so that z axis is actual z axis and cell # is particle #
imsz=size(iStackTmp{1},1);
IMz = zeros(imsz,imsz,nFr);
for ii = 1:nPsf
    imPsfStack{ii} = IMz;%preallocate
    for jj = 1:nFr
        imPsfStack{ii}(:,:,jj) = iStackTmp{jj}(:,:,ii);
    end
end
%--------------------------------------------
function [drift] = getImDrift(templateIm,featureIm,windowSize)
% spatial cross correlation function = SCCF
if ~exist('windowSize','var')
  windowSize= 10; %seems like a pretty reasonable number
end

%calculate the correlation function
% C is the image of the correlation function. 
%zeroCoord is the [i,j] coordinate corresponding to zero displacement in
% the correlation function 
[C,zeroCoord] = corrfunc(templateIm,featureIm);


% careful with (i,j) vs (x,y)!!
[corrMaxPos corAmplitude] =   getPeakPosCentroid(C,zeroCoord,windowSize); 
drift = zeroCoord  - corrMaxPos;
%-------------------------------------------------------------
function [G,zeroCoord] = corrfunc(template,feature)
% July 9, 2003
% David Kolin
% 18/03/11
% Seamus Holden
% 1)Minor modification 2011 SH to output zeroCoordinate in ICS image
% 2) 110318 Now a very heavily modified version of original function. Does cross not auto correlation
% NB:template and feature should be same size
% zeroCoord is (x,y) coordinates, not (i,j)!
template=double(template);
feature=double(feature); 

% Calculates 2D correlation functions for each z slice of a given 3D matrix (imgser)
% Output as a 3D matrix with same dimensions as imgser
G = zeros(size(template)); % Preallocates matrix

% Calculates corr func
%autocorrelation:
%%G = (fftshift(real(ifft2(fft2(template).*conj(fft2(template)))))) ...
%%        /(mean(template(:))^2*size(template,1)*size(template,2) ) ...
%%      -1;

%cross correlation
G = (fftshift(real(ifft2(fft2(template).*conj(fft2(feature))))))/...
		( (mean(template(:))*mean(feature(:))) * size(template,1)*size(template,2) ) ...
		- 1;
% SH mod
% make sure that you know where the zero coordinate is from the DFT
% so that we can calculate absolute drift
imsize = size(template);
zeroCoordX = (floor(imsize(2)/2)+1);
zeroCoordY = (floor(imsize(1)/2)+1);
zeroCoord = [ zeroCoordX,zeroCoordY];

%----------------------------------------------
function psfCtr = getPsfCtr(imPsf,windowSz)

imsz = size(imPsf,1);
imctrX= floor(imsz/2)+1;
imctr=[imctrX,imctrX];

nFr = size(imPsf,3);
for ii=1:nFr
    [psfCtr(ii,:)]= getPeakPosCentroid(imPsf(:,:,ii), imctr, windowSz);
end
%---------------------------------------------------
function posListCorr = crossCorrPsfPos(imPsfInFocus,posListCtrd,cropRatio,windowSize)
%cross correlate in focus PSFs to 1st psf to find updated centre pos

%windowSize = round(HLFBOX*cropRatio/2);%for the centroid finding in th CCF image. Ie PSF to PSF shift can be big but not huge
im1 = imPsfInFocus(:,:,1);
nMol= size(posListCtrd,1);
for ii= 1:nMol
    imN = imPsfInFocus(:,:,ii);
    [psfShift(ii,:)] = getImDrift(im1,imN, windowSize);
    posListCorr(ii,:) = posListCtrd(ii,:)+psfShift(ii,:)/cropRatio;
end



%---------------------------------------------------------
function imPsfCtrd= recentre(imPsf,psfCtr,pixSzOut, MAXSHIFT);
%shift each frame so CoM is at middle of image (ie floor(sizeX/2)+1)

nFr = size(psfCtr,1);
szY = size(imPsf,1);
ctrY = floor(szY/2)+1;
%d=sqrt(sum((psfCtr-imCtr).^2,2));
%dNm=d*pixSzOut;
imCtr=[ctrY,ctrY];
d = psfCtr-imCtr;
biggestShift = ceil(max(abs(d)));

if biggestShift*pixSzOut>=MAXSHIFT
    error('Shift in PSF centroid position too large for realignment. Reduce Z-range, increase MAXSHIFT, or increase quality of data.');
else
    cropLim = floor(szY/2-biggestShift);
    x = round(psfCtr(1));
    y = round(psfCtr(2));
    imPsfCtrd= imPsf(y-cropLim:y+cropLim, x-cropLim:x+cropLim,:);
end

%---------------------------------------------------------
function imstackAvg=zFrameAverage(imstack,framePerZPos)

imstack = double(imstack);

sz = size(imstack);
nFr = sz(3);
if mod(nFr,framePerZPos)~=0
    error('Total number of frames not divisible by number per Z-position');
else
    %average the frames
    nFr2 = nFr/framePerZPos;
    imstackAvg = zeros(sz(1),sz(2),nFr2);
    imAvg = zeros(sz(1),sz(2));
    for ii = 1:nFr2
        imAvg=imAvg*0;
        for jj=1:framePerZPos
            curFrame = (ii-1)*framePerZPos+jj;
            imAvg = imAvg+imstack(:,:,curFrame)/framePerZPos;
        end
        imstackAvg(:,:,ii)=imAvg;
    end
end
%---------------------------------------------------------
function imstack2= interpZ(imstack, interpZFactor);
%assumes imstack of type double

%interpolated z
nFr = size(imstack,3);
nY= size(imstack,1);
nX= size(imstack,2);

z1 = 1:nFr;
z2 = 1:(1/interpZFactor):nFr;
nFr2 = numel(z2);
imstack2 = zeros(nY,nX,nFr2);
%interpolate each pixel
for ii = 1:nY
    for jj = 1:nX
       i1 = squeeze(imstack(ii,jj,:));
       i2 = interp1(z1,i1,z2,'pchip');
       imstack2(ii,jj,:) = reshape(i2,1,1,nFr2);
    end
end