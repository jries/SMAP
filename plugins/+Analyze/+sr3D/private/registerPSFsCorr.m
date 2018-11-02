function [outstack,shiftedim]=registerPSFsCorr(imstackin,posfocus)
% imstack=imstack(:,:,:,1:2);
s=size(imstackin);
imstack=zeros(s(1),s(2),s(3),s(4)+1);
imstack(:,:,:,2:end)=imstackin;
imstack(:,:,:,1)=sum(imstackin,4);
posListCtrd=zeros(size(imstack,4),2)+round((s(1)+1)/2);
% posfocus=12;
imPsfInFocus=squeeze(imstack(:,:,posfocus,:));
cropRatio=1;
windowSize=10;
posListCorr = crossCorrPsfPos(imPsfInFocus,posListCtrd,cropRatio,windowSize);
shift=posListCorr-round((s(1)+1)/2);
meanim=zeros(s(1:3));
shiftedim=zeros(size(imstackin));
for k=1:s(4)
    shiftedim(:,:,:,k)=imtranslate(imstack(:,:,:,k+1),-shift(k+1,[1 2]));
    meanim=meanim+shiftedim(:,:,:,k);
end

% posListCorr = crossCorrPsfPos(squeeze(shiftedim(:,:,posfocus,:)),posListCtrd(1:end-1,:),cropRatio,windowSize)
outstack=meanim/s(4);

 imagesc([squeeze(sum(imstackin(:,:,posfocus,:),4)) meanim(:,:,posfocus)])
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



function [pos, amplitude] = getCentroid(im);
im = im-min(im(:));%offset the image to avoid baseline "zero" values contributing to the avg

[sizey sizex] = size(im);
totalIntensity = sum(im(:));
posx = sum(sum(im,1).*[1:sizex]) ./ totalIntensity;
posy = sum(sum(im,2).*[1:sizey]') ./ totalIntensity;

pos = [posx,posy];
amplitude = im(round(posy),round(posx));

