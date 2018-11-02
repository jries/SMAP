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

function f_scatter3(totalPSFfits,useFidCorrections)
%% Ask for plotting parameters
scrsz = get(0,'ScreenSize');
% c_map = hot(256);
useFidCorrections = logical(useFidCorrections);

dlg_title = 'Input plotting parameters';
prompt = {'Immersion lens index of refraction', ...
    'Sample media index of refraction', ...
    'Size of scatter points',...
    'Temporal color map (0=off, other=on)',...
    'z-range to plot (nm)',...
    'Frame range to plot'};    %,...
def = {'1.518','1.33','10','0','[-2000 2000]',...
    mat2str([min(totalPSFfits(:,1)) max(totalPSFfits(:,1))])};    %, ... 
inputdialog = inputdlg(prompt,dlg_title,1,def);

% useFidCorrections = str2double(inputdialog{1});
% powerAtObjective = str2double(inputdialog{5})/1000;


% nmPerPixel = 125.78;    % was 160 for 8b back
% scaleBarLength = 1000;  % nm
% pixelSize = 2;          % size of pixels in reconstructed image in nm
% border = 500;           % plot extra region around the cells (size of extra region in nm)
% lambda = 615;           % nm, was 527
% NA = 1.4;               % numerical aperture
nOil = str2double(inputdialog{1});        % index of immersion oil
nSample = str2double(inputdialog{2});     % index of refraction of sample
scatterSize = str2double(inputdialog{3});
useTimeColors = logical(str2double(inputdialog{4}));
zRange = str2num(inputdialog{5});
frameRange = str2num(inputdialog{6});


%% Now re-evaluate the goodness of the fits

% Conditions for fits (play with these):
% (1) Amplitude of both lobes > 0
% (2) All locations x1,y1, x2,y2 lie inside area of small box
% (3) All sigmas need to be > sigmaBound(1) and < sigmaBound(2)
% (4) Distance between lobes needs to be > lobeDist(1) pixels and < lobeDist(2) pixels
% (5) Make sure amplitudes are within 100% of one another
% (6) Make sure totalFitError/(total number of photons) is within the fitErrorRange

% fitErrorCol = 16;
% goodFitFlagCol = 17;
% numPhotonCol = 21;
% lobeDistCol = 22;
% ampRatioCol = 23;
% sigmaRatioCol = 24;


% for i = 1:size(totalPSFfits,1)
% 
%     if totalPSFfits(i,goodFitFlagCol) == -1001 || ...
%             totalPSFfits(i,goodFitFlagCol) == -1002 || ...
%             totalPSFfits(i,goodFitFlagCol) == -1003
% 
%         continue
% 
%     elseif totalPSFfits(i,sigmaRatioCol) > sigmaRatioLimit;
% 
%         totalPSFfits(i,goodFitFlagCol) = -1004;
% 
%     elseif totalPSFfits(i,lobeDistCol) < lobeDistBounds(1) || totalPSFfits(i,lobeDistCol) > lobeDistBounds(2)
% 
%         totalPSFfits(i,goodFitFlagCol) = -1005;
% 
%     elseif totalPSFfits(i,ampRatioCol) > ampRatioLimit;
% 
%         totalPSFfits(i,goodFitFlagCol) = -1006;
% 
%     elseif totalPSFfits(i,fitErrorCol)*conversionFactor/totalPSFfits(i,numPhotonCol) > fitErrorRange(2)  || ...
%             totalPSFfits(i,fitErrorCol)*conversionFactor/totalPSFfits(i,numPhotonCol) < fitErrorRange(1)
% 
%         totalPSFfits(i,goodFitFlagCol) = -1007;
%     else
%         totalPSFfits(i,goodFitFlagCol) = 3;
%     end
% 
% end

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
frameNum = frameNum(validPoints);

%% plot 3D scatterplot of localizations with white light

f = figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','k','renderer','opengl');

if useTimeColors == 0
    plot3(xLoc,yLoc,zLoc,'.','MarkerSize',scatterSize,...
        'Color',[1 1 0]);
%     scatter3(xLoc,yLoc,zLoc,scatterSize,[1 1 0],'filled');
else
    %         scatter3(xLoc(a),yLoc(a),zLoc(a),scatterSize,frameNum(1):frameNum(length(frameNum)),'filled')
    markerColors = jet(frameNum(length(frameNum))-frameNum(1)+1);
    %     for a = 1:length(frameNum)
    %         scatter3(xLoc(a)-min(xLoc(:)),yLoc(a)-min(yLoc(:)),zLoc(a),scatterSize,'filled',...
    %             'MarkerFaceColor', markerColors(frameNum(a)-frameNum(1)+1,:),...
    %             'MarkerEdgeColor', markerColors(frameNum(a)-frameNum(1)+1,:));
    %     end

    for a = 1:length(frameNum)
        plot3(xLoc(a),yLoc(a),zLoc(a),'.','MarkerSize',scatterSize,...
            'Color',markerColors(frameNum(a)-frameNum(1)+1,:));
        hold on;
%         scatter3(xLoc(a),yLoc(a),zLoc(a),scatterSize,'filled',...
%             'MarkerFaceColor', markerColors(frameNum(a)-frameNum(1)+1,:),...
%             'MarkerEdgeColor', markerColors(frameNum(a)-frameNum(1)+1,:));
    end
end
grid on;
axis image vis3d;
%     xlim([min(xLoc(:)) max(xLoc(:))]);
%     ylim([min(yLoc(:)) max(yLoc(:))]);
%     xlim([min(xLoc(:)) max(xLoc(:))]-min(xLoc(:)));
%     ylim([min(yLoc(:)) max(yLoc(:))]-min(yLoc(:)));
% xlim([min(xLoc(:)) max(xLoc(:))]);
% ylim([min(yLoc(:)) max(yLoc(:))]);
xlabel('x (\mum)');ylabel('y (\mum)');zlabel('z (\mum)');

title([num2str(length(xLoc)) ' localizations'],'color','w');
set(gca,'color','k');
set(gca,'xcolor','w');set(gca,'ycolor','w');set(gca,'zcolor','w');


end
