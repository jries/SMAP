% Copyright (c)2013, The Board of Trustees of The Leland Stanford Junior
% University. All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
% Redistributions of source code must retain the above copyright notice, 
% this list of conditions and the following disclaimer.
% Redistributions in binary form must reproduce the above copyright notice, 
% this list of conditions and the following disclaimer in the documentation 
% and/or other materials provided with the distribution.
% Neither the name of the Leland Stanford Junior University nor the names 
% of its contributors may be used to endorse or promote products derived 
% from this software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS 
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% plot 2D images of molecule locations

function f_hist2(totalPSFfits,useFidCorrections)
%% Ask for plotting parameters
scrsz = get(0,'ScreenSize');
% c_map = hot(256);
useFidCorrections = logical(useFidCorrections);

dlg_title = 'Input plotting parameters';
prompt = {'Immersion lens index of refraction', ...
    'Sample media index of refraction', ...
    'Histogram bin size (nm)',...
    'z-range to plot (nm)',...
    'Frame range to plot'};    %,...
def = {'1.518','1.33','25','[-2000 2000]',...
    mat2str([min(totalPSFfits(:,1)) max(totalPSFfits(:,1))])};    %, ... 
inputdialog = inputdlg(prompt,dlg_title,1,def);

% useFidCorrections = str2double(inputdialog{1});
% powerAtObjective = str2double(inputdialog{5})/1000;


% nmPerPixel = 125.78;    % was 160 for 8b back
% scaleBarLength = 1000;  % nm
% binSize = 2;          % size of pixels in reconstructed image in nm
% border = 500;           % plot extra region around the cells (size of extra region in nm)
% lambda = 615;           % nm, was 527
% NA = 1.4;               % numerical aperture
nOil = str2double(inputdialog{1});        % index of immersion oil
nSample = str2double(inputdialog{2});     % index of refraction of sample
binSize = str2double(inputdialog{3});
zRange = str2num(inputdialog{4});
frameRange = str2num(inputdialog{5});


%% load valid xyz locations

%     goodFits = false(size(totalPSFfits,1),1);

goodFits = totalPSFfits(:,17) > 0; % totalPSFfits(:,11) > -inf;
goodFits = goodFits & totalPSFfits(:,1) >= frameRange(1) & totalPSFfits(:,1) <= frameRange(2);
% goodFits = goodFits & totalPSFfits(:,numPhotonCol) >= numPhotonRange(1) & totalPSFfits(:,numPhotonCol) <= numPhotonRange(2);

if useFidCorrections
    goodFits = goodFits & totalPSFfits(:,30) >= zRange(1) & totalPSFfits(:,30) <= zRange(2);
    xLoc = totalPSFfits(goodFits,28)/1000;  % convert from nm to um
    yLoc = totalPSFfits(goodFits,29)/1000;
    zLoc = totalPSFfits(goodFits,30)/1000;

else
    goodFits = goodFits & totalPSFfits(:,27) >= zRange(1) & totalPSFfits(:,27) <= zRange(2);
    xLoc = totalPSFfits(goodFits,25)/1000;
    yLoc = totalPSFfits(goodFits,26)/1000;
    zLoc = totalPSFfits(goodFits,27)/1000;

end

binSize = binSize/1000;
% simple correction for aberrations if sample is near coverslip
zLoc = zLoc * nSample/nOil;
numPhotons = totalPSFfits(goodFits,21);
%     meanBkgnd = totalPSFfits(goodFits,15)*conversionFactor;
% meanBkgnd = totalPSFfits(goodFits,15);  % output from template match is already in units of photons.
frameNum = totalPSFfits(goodFits,1);

%% ask user what region to plot in superresolution image

h=figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');

% plot is faster than scatter
plot(xLoc,yLoc,'.','MarkerSize',1);
% xlim([min(xLoc(:)) max(xLoc(:))]);
% ylim([min(yLoc(:)) max(yLoc(:))]);
xlabel('x (\mum)');ylabel('y (\mum)');
axis image ij;

ROI = imrect(gca,[min(xLoc(:)) min(yLoc(:)) max(xLoc(:))-min(xLoc(:)) max(yLoc(:))-min(yLoc(:))]);

title({'Double-click to choose region that will be plotted in 3D scatterplot' ...
    mat2str(ROI.getPosition)});
addNewPositionCallback(ROI,@(p) title({'Double-click to choose region that will be plotted in 3D scatterplot' ...
    ['[xmin ymin width height] = ' mat2str(p,3)]}));
% make sure rectangle stays within image bounds
fcn = makeConstrainToRectFcn('imrect',get(gca,'XLim'),get(gca,'YLim'));
setPositionConstraintFcn(ROI,fcn);
ROI = wait(ROI);
clear avgImg fcn
close(h);

%% filter out localizations outside of ROI

validPoints = xLoc>=ROI(1) & xLoc<=ROI(1)+ROI(3) & yLoc>ROI(2) & yLoc<=ROI(2)+ROI(4) & numPhotons>0;

xLoc = xLoc(validPoints);
yLoc = yLoc(validPoints);
zLoc = zLoc(validPoints);
% meanBkgnd = meanBkgnd(validPoints);
% frameNum = frameNum(validPoints);

%% histogram localizations
xEdges = min(xLoc):binSize:(max(xLoc)+binSize);
yEdges = min(yLoc):binSize:(max(yLoc)+binSize);
xPx = floor((xLoc-min(xLoc))/binSize)+1;
yPx = floor((yLoc-min(yLoc))/binSize)+1;

zHist = zeros(length(yEdges)-1,length(xEdges)-1);
brightHist = zeros(length(yEdges)-1,length(xEdges)-1);
pointIdx = sub2ind(size(zHist),yPx,xPx);
%zRange = [min(zLoc(zLoc>=-1000)) max(zLoc(zLoc<=1000))];
for idx=unique(pointIdx)'
    points = idx == pointIdx;
    zHist(idx) = median(zLoc(points));
    brightHist(idx) = length(zLoc(points));
end
zHist(brightHist == 0) = NaN;


%% create histogram colorbar
brightMax = round(prctile(brightHist(:),[0.5 99.5]));
brightMin = brightMax(1);
brightMax = brightMax(2);
% brightMax = max(brightHist(:));     %(mean(brightHist(:))+2*std(brightHist(:)));
brightHistImg = (brightHist-brightMin)/(brightMax-brightMin);
brightHistImg(brightHistImg>1) = 1;
brightHistImg(brightHistImg<0) = 0;

zRange = prctile(zLoc,[2 98]);  %[min(zLoc) max(zLoc)];
brightBar = 1:max(brightHist(:));
z_colorBar = linspace(zRange(1),zRange(2),length(yEdges)-1);
[brightBar, zBar] = meshgrid(brightBar,z_colorBar);
histBar = ind2rgb(round((zBar-zRange(1))/(zRange(2)-zRange(1))*255)+1,jet(256));
histBar = rgb2hsv(histBar);

brightBar = (brightBar-brightMin)/(brightMax-brightMin);
brightBar(brightBar>1) = 1;
brightBar(brightBar<0) = 0;
histBar(:,:,3) = histBar(:,:,3) .* brightBar;
%histBar(:,:,3) = brightBar/max(brightBar(:));
histBar = hsv2rgb(histBar);

%% generate histogram image
histImg = ind2rgb(round((zHist-zRange(1))/(zRange(2)-zRange(1))*255)+1,jet(256));
histImg = rgb2hsv(histImg);
histImg(:,:,3) = histImg(:,:,3) .* brightHistImg;
%histImg(:,:,3) = brightHist/max(brightHist(:));
histImg = hsv2rgb(histImg);



%% plot 2D histograms

figure('Position',[(scrsz(3)-1280)/2+1 (scrsz(4)-720)/2 1280 720],'color','w');
% set(gcf,'DefaultTextFontSize',12,'DefaultAxesFontSize',12);
imagesc(xEdges,yEdges,zHist,zRange);axis image; 
colormap([0 0 0; jet(63)]);
colorbar;
title(['Histogram of median z-position: ' num2str(length(zLoc)) ' localizations']);
xlabel('x (\mum)');
ylabel('y (\mum)');

% figure('Position',[(scrsz(3)-1280)/2+1 (scrsz(4)-720)/2 1280 720],'color','w');
% set(gcf,'DefaultTextFontSize',12,'DefaultAxesFontSize',12);
% imagesc(xEdges,yEdges,brightHist);axis image; colormap hot; colorbar;
% title(['Localization density: ' num2str(length(zLoc)) ' localizations'])
% xlabel('x (\mum)');
% ylabel('y (\mum)');

figure('Position',[(scrsz(3)-1280)/2+1 (scrsz(4)-720)/2 1280 720],'color','w');
% set(gcf,'DefaultTextFontSize',12,'DefaultAxesFontSize',12);
h=subplot(1,8,1:7);
image(xEdges,yEdges,histImg);axis image;
title(['Histogram of number of localizations and median z-position: ' num2str(length(zLoc)) ' localizations']);
xlabel(h,'x (\mum)');
ylabel(h,'y (\mum)');

h=subplot(1,8,8);
image(1:max(brightHist(:)),zRange,histBar);axis xy;
xlabel(h,'# localizations');
ylabel(h,'z (\mum)');
set(gca,'YAxisLocation','right');
