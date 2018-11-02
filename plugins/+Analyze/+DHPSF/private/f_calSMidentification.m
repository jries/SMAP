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

% function [templateFrames, ROI, dataFile, dataPath, darkFile, logFile,...
%             logPath, EMGain, templateLocs] = f_calSMidentification(calFile,calBeadIdx,...
%             p.templateFile, p.boxRadius,channel,p.sigmaBounds,p.gaussianFilterSigma,p.minDistBetweenSMs)
function [templateFrames, ROI, dataFile, dataPath, templateLocs,outputFilePrefix,imall] = f_calSMidentification(obj,p)
        
        
ploton=false;
% f_calSMidentification is a module in easy_dhpsf that prepares the
% templates from f_calDHPSF and uses them to generate a series of template
% matches. These are then used to judge an appropriate threshold for 
% f_fitSMs. This module also sets the file and other parameters used for 
% f_fitSMs, as well as some parameters for f_trackFiducials.

% Instrument Specific Parameters
channel='0';
calBeadIdx=1; %XXX extend
% dlg_title = 'Set EM Gain';
% prompt = {  'EM Gain (1 if no gain):' };
% def = {'300'};
% 
% num_lines = 1;
% inputdialog = inputdlg(prompt,dlg_title,num_lines,def);
% 
% if isempty(inputdialog)
%     error('User cancelled the program')
% end
% 
% EMGain = str2double(inputdialog{1});
% if EMGain < 1 || isnan(EMGain)
%     warning('EMGain should be >= 1. Setting to 1...');
%     EMGain = 1;
% end

frameNum = 1;
scrsz = get(0,'ScreenSize');
% Options for lsqnonlin
options = optimset('FunValCheck','on','Diagnostics','off','Jacobian','on', 'Display', 'off');
%    'FinDiffType','central','DerivativeCheck','on');

%% ask user for relevant datafiles

% [dataFile, dataPath] = uigetfile({'*.tif';'*.*'},...
%     'Open SMACM image stack(s) for data processing',...
%     'MultiSelect', 'on');
% if isequal(dataFile,0)
%     error('User cancelled the program');
% end
% 
% if ischar(dataFile)
%     dataFile = cellstr(dataFile);
% end
[dataPath,fi,ex]=fileparts(p.dataFile);
dataPath=[dataPath filesep];
dataFile={[fi ex]};

for stack = 1:length(dataFile)

    fileInfo = imfinfo([dataPath dataFile{stack}]);
    numFrames = length(fileInfo);
    frames = 1:numFrames;
    imgHeight = fileInfo(1).Height;
    imgWidth = fileInfo(1).Width;
    
    if stack == 1
        
%         [p.templateFile, templatePath] = uigetfile({'*.tif;*.mat';'*.*'},'Open image stack or MATLAB *.mat file of a DH-PSF template');
%         if isequal(p.templateFile,0)
%             error('User cancelled the program');
%         end
        
        
        if strcmp(p.templateFile(length(p.templateFile)-2:length(p.templateFile)),'tif')
            
            templateInfo = imfinfo(p.templateFile);
            if templateInfo(1).Height ~= templateInfo(1).Width
                error('Template is not square');
            end
            templateSize = templateInfo(1).Height;
        else
            load(p.templateFile);
            templateSize = size(template,2);
        end
        
        load(p.calFile);
        
        % compute templates to use based upon fitted DG angles
        templateFrames = interp1(squeeze(meanAngles(1,calBeadIdx,goodFit_forward)),...
            1:length(meanAngles(1,calBeadIdx,goodFit_forward)),-60:30:90,'nearest'); %90:-30:-60
 
        % if some angles are out of range, use the ends of the template
        % stack. first determine whether frames are increasing or
        % decreasing so the the 'ends' are chosen correctly.
        
        if any(isnan(templateFrames))
            if nanmean(diff(templateFrames)) >= 0
                endFrames = {1 sum(goodFit_forward)};
                df=-1;
            elseif nanmean(diff(templateFrames)) < 0
                endFrames = {sum(goodFit_forward) 1};
                df=1;
            else
                endFrames = nan;
            end
            if isnan(templateFrames(1))
                templateFrames(1) = endFrames{1};
            end
            if isnan(templateFrames(end))
                templateFrames(end) = endFrames{2};
            end
        end

        for m=length(templateFrames)-1:-1:2
            if isnan(templateFrames(m))
                templateFrames(m)=(templateFrames(m+1)+templateFrames(m-1))/2;
            end
                
        end
%         temp = inputdlg({['Input sets of frames corresponding to each template or the index of templates to use (ex. ' ...
%             mat2str(templateFrames) ')']},...
%             'Input template numbers',1, ...
%             {mat2str(templateFrames)}); ...     % {'[8:4:32]'});
%         templateFrames = str2num(temp{1});
        
%         logPath = 0;
%         logFile = 0;
%         [logFile, logPath] = uigetfile({'*.dat';'*.*'},...
%             'Open sequence log file(s) corresponding to image stack(s) (optional: hit cancel to skip)',...
%             'MultiSelect', 'on');
%         if isequal(logPath,0)
%             logFile = 'not specified';
%         end
%         if ischar(logFile)
%             logFile = cellstr(logFile);
%         end

%         [darkFile, darkPath] = uigetfile({'*.tif';'*.*'},'Open image stack with dark counts (same parameters as SMACM data)');
%         if isequal(darkFile,0)
%             error('User cancelled the program');
%         end
        
    end
    
    %% create output log filenames
    if channel == '0'
        outputFilePrefix{stack} = [dataPath dataFile{stack}(1:length(dataFile{stack})-4) filesep 'threshold ' ...
            datestr(now,'yyyymmdd HHMM') filesep];
        mkdir(outputFilePrefix{stack});
    elseif channel == 'G'
        outputFilePrefix{stack} = [dataPath dataFile{stack}(1:length(dataFile{stack})-4)  filesep 'reflected' filesep 'threshold ' ...
            datestr(now,'yyyymmdd HHMM') filesep];
        mkdir(outputFilePrefix{stack});
    elseif channel == 'R'
        outputFilePrefix{stack} = [dataPath dataFile{stack}(1:length(dataFile{stack})-4) filesep 'transmitted' filesep 'threshold ' ...
            datestr(now,'yyyymmdd HHMM') filesep];
        mkdir(outputFilePrefix{stack});
    end
    
    if stack==1
        %% Compute darkAvg counts
        
%         if ~isequal(darkFile,0)
%             darkFile = [darkPath darkFile];
%             % Computes average of darkAvg frames for background subtraction
%             darkFileInfo = imfinfo(darkFile);
%             numDarkFrames = length(darkFileInfo);
%             darkAvg = zeros(darkFileInfo(1).Height,darkFileInfo(1).Width);
%             for frame = 1:numDarkFrames
%                 darkAvg = darkAvg + double(imread(darkFile,frame,'Info',darkFileInfo));
%             end
%             darkAvg = darkAvg/numDarkFrames;
%             if ~isequal(size(darkAvg),[imgHeight imgWidth])
%                 warning('Dark count image and data image stack are not the same size. Resizing darkAvg count image...');
%                 darkAvg = imresize(darkAvg,[imgHeight imgWidth]);
%             end
%         else
%             darkAvg = 0;
%         end
%         clear darkFileInfo;
        darkAvg=p.offsetADU;
        
 
        
        %% opens and processes templates of DH-PSF for template matching
        
        if strcmp(p.templateFile(length(p.templateFile)-2:length(p.templateFile)),'tif')
            
            numTemplates = size(templateFrames,1);
            templateColors = jet(numTemplates);
            template = zeros(numTemplates,templateSize,templateSize);
            templateLocs = zeros(numTemplates,5);
            fitParam = zeros(1,8);
            [xIdx, yIdx] = meshgrid(1:templateSize,1:templateSize);
            
%             if ploton
                hTemplate=figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
                set(hTemplate,'Visible','off');
%             end
            for a=1:numTemplates
                for b=templateFrames(a,1):templateFrames(a,2)
                    template(a,:,:) = squeeze(template(a,:,:)) + ...
                        double(imread([templatePath p.templateFile],b,'Info',templateInfo));
                end
                % make minimum count level in template equal to 0
                template(a,:,:) = template(a,:,:) - min(min(template(a,:,:)));
                % normalize energy contained (sum of all counts) in the template
                template(a,:,:) = template(a,:,:) / sum(sum(template(a,:,:)));
                % finally, make mean of template equal to 0
                template(a,:,:) = template(a,:,:) - mean(mean(template(a,:,:)));
                
                % find two largest peaks in template
                [tempY, tempX] = ind2sub([templateSize templateSize],find(imregionalmax(template(a,:,:))));
                temp = sortrows([tempX tempY template(sub2ind([numTemplates templateSize templateSize],a*ones(length(tempX),1),tempY,tempX))],-3);
                
                % [amp1 amp2 xMean1 yMean1 xMean2 yMean2 sigma1 sigma2]
                fitParam(3) = temp(1,1);
                fitParam(4) = temp(1,2);
                fitParam(5) = temp(2,1);
                fitParam(6) = temp(2,2);
                fitParam(1) = temp(1,3);
                fitParam(2) = temp(2,3);
                fitParam(7) = 1.8;
                fitParam(8) = 1.8;
                lowerBound = [0 0 1 1 1 1 p.sigmaBounds(1) p.sigmaBounds(1)];
                upperBound = [max(max(template(a,:,:))) max(max(template(a,:,:))) ...
                    templateSize templateSize templateSize templateSize ...
                    p.sigmaBounds(2) p.sigmaBounds(2)];
                
                % Fit with lsqnonlin
                fitParam = lsqnonlin(@(x) ...
                    f_doubleGaussianVector(x,squeeze(template(a,:,:)),0,xIdx,yIdx),...
                    fitParam,lowerBound,upperBound,options);
                
                templateLocs(a,1:2) = fitParam(3:4);
                templateLocs(a,3:4) = fitParam(5:6);
                % calculate rough angle between peaks
                templateLocs(a,5) = 180/pi*atan2(templateLocs(a,2)-templateLocs(a,4), ...
                    templateLocs(a,1)-templateLocs(a,3));
                
%                 if ploton
                    subplot(1,numTemplates,a);imagesc(squeeze(template(a,:,:)));
                    axis image;colormap hot;colorbar;
                    hold on;
                    plot(templateLocs(a,1),templateLocs(a,2),'.','MarkerEdgeColor', templateColors(a,:));
                    plot(templateLocs(a,3),templateLocs(a,4),'.','MarkerEdgeColor', templateColors(a,:));
                    title({['Template ' mat2str(templateFrames(a,:))] ...
                        ['Angle = ' num2str(templateLocs(a,5)) ' deg']});
%                 end
            end
            imwrite(frame2im(getframe(hTemplate)),[outputFilePrefix{stack} 'templates.tif']);
            clear templateInfo tempX tempY temp xIdx yIdx;

        else
            templateFrames = round( templateFrames');
            numTemplates = size(templateFrames,1);
            templateColors = jet(numTemplates);
            templateLocs = zeros(numTemplates,5);
            fitParam = zeros(1,8);
            [xIdx, yIdx] = meshgrid(1:templateSize,1:templateSize);
            
            
            hTemplate=figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
            for a=1:numTemplates
                
                % make minimum count level in template equal to 0
                template(templateFrames(a),:,:) = template(templateFrames(a),:,:)...
                    - min(min(template(templateFrames(a),:,:)));
                % normalize energy contained (sum of all counts) in the template
                template(templateFrames(a),:,:) = template(templateFrames(a),:,:)...
                    / sum(sum(template(templateFrames(a),:,:)));
                % finally, make mean of template equal to 0
                template(templateFrames(a),:,:) = template(templateFrames(a),:,:)...
                    - mean(mean(template(templateFrames(a),:,:)));
                
                
                % find two largest peaks in template
                [tempY, tempX] = ind2sub([templateSize templateSize],find(imregionalmax(template(templateFrames(a),:,:))));
                temp = sortrows([tempX tempY template(sub2ind(size(template),...
                    templateFrames(a)*ones(length(tempX),1),tempY,tempX))],-3);
                
                % [amp1 amp2 xMean1 yMean1 xMean2 yMean2 sigma1 sigma2]
                fitParam(3) = temp(1,1);
                fitParam(4) = temp(1,2);
                fitParam(5) = temp(2,1);
                fitParam(6) = temp(2,2);
                fitParam(1) = temp(1,3);
                fitParam(2) = temp(2,3);
                fitParam(7) = mean(p.sigmaBounds);
                fitParam(8) = mean(p.sigmaBounds);
                lowerBound = [0 0 1 1 1 1 p.sigmaBounds(1) p.sigmaBounds(1)];
                upperBound = [max(max(template(templateFrames(a),:,:))) max(max(template(templateFrames(a),:,:))) ...
                    templateSize templateSize templateSize templateSize ...
                    p.sigmaBounds(2) p.sigmaBounds(2)];
                
                % Fit with lsqnonlin
                fitParam = lsqnonlin(@(x) ...
                    f_doubleGaussianVector(x,squeeze(template(templateFrames(a),:,:)),0,xIdx,yIdx),...
                    fitParam,lowerBound,upperBound,options);
                
                templateLocs(a,1:2) = fitParam(3:4);
                templateLocs(a,3:4) = fitParam(5:6);
                % calculate rough angle between peaks
                templateLocs(a,5) = 180/pi*atan2(templateLocs(a,2)-templateLocs(a,4), ...
                    templateLocs(a,1)-templateLocs(a,3));
                
                subplot(1,numTemplates,a);
                imagesc(squeeze(template(templateFrames(a),:,:)));
                axis image;colormap hot;colorbar;
                hold on;
                plot(templateLocs(a,1),templateLocs(a,2),'.','MarkerEdgeColor', templateColors(a,:));
                plot(templateLocs(a,3),templateLocs(a,4),'.','MarkerEdgeColor', templateColors(a,:));
                title({['Frames ' mat2str(templateFrames(a,:))] ...
                    ['Angle = ' num2str(templateLocs(a,5)) ' deg']});
                %         drawnow
                %         pause(1)
            end
            imwrite(frame2im(getframe(hTemplate)),[outputFilePrefix{stack} 'templates.tif']);
            clear templateInfo tempX tempY temp xIdx yIdx;

        end
        close(hTemplate); % closes template figure
        %% user picks ROI
        % pick region of interest by reading first frame and having user select
        % region
        
%         % Compute average image
%         avgImg = zeros(imgHeight,imgWidth);
%         for a = 1:min(numFrames,200)
%             avgImg = avgImg + double(imread([dataPath dataFile{stack}],frames(a),'Info',fileInfo)) - darkAvg;
%         end
%         avgImg = avgImg/200;
%         
%         hROI = figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
%         imagesc(avgImg);axis image;colormap hot;
%         if channel == 'G'
%             ROI = imrect(gca,[1 1 270 270]);
%         elseif channel == 'R'
%             ROI = imrect(gca,[243 243 270 270]);
%         else
%             ROI = imrect(gca,[1 1 size(avgImg,1) size(avgImg,2)]);
%         end
%         
%         % ROI = imrect(gca,[1 1 128 128]);
%         title({'Shape box and double-click to choose region of interest for PSF extraction' ...
%             ['[xmin ymin width height] = ' mat2str(ROI.getPosition)]...
%             'The displayed image is the average of the first 200 frames'});
%         addNewPositionCallback(ROI,@(p) title({'Shape box and double-click to choose region of interest for PSF extraction' ...
%             ['[xmin ymin width height] = ' mat2str(p,3)]...
%             'The displayed image is the average of the first 200 frames'}));
%         % make sure rectangle stays within image bounds
%         fcn = makeConstrainToRectFcn('imrect',get(gca,'XLim'),get(gca,'YLim'));
%         setPositionConstraintFcn(ROI,fcn);
%         ROI = round(wait(ROI));
%         % make sure ROI is an even number of pixels
%         close(hROI) % closes ROI selection
ROI=[1 1 imgHeight imgWidth];
        if mod(ROI(3),2)==1
            ROI(3) = ROI(3)-1;
        end
        if mod(ROI(4),2)==1
            ROI(4) = ROI(4)-1;
        end
%         %ROI = [84 127 128 130];
        cropWidth = ROI(3);
        cropHeight = ROI(4);

        
        %% prepare template for template matching
        
        % pad template to same size as input
        templatePad = zeros(numTemplates,cropHeight,cropWidth);
        templateFT = zeros(numTemplates,cropHeight,cropWidth);
        for a=1:numTemplates
            
            if strcmp(p.templateFile(length(p.templateFile)-2:length(p.templateFile)),'tif')
                templatePad(a,:,:) = padarray(squeeze(template(a,:,:)),...
                    [(cropHeight-size(template,2))/2 ...
                    (cropWidth-size(template,3))/2],min(min(template(a,:,:))));
            else
                templatePad(a,:,:) = padarray(squeeze(template(templateFrames(a),:,:)),...
                    [(cropHeight-size(template,2))/2 ...
                    (cropWidth-size(template,3))/2],min(min(template(templateFrames(a),:,:))));
            end
            
            % multiplying by conjugate of template in FT domain is squivalent
            % to flipping the template in the real domain
            templateFT(a,:,:) = conj(fft2(squeeze(templatePad(a,:,:))));
        end
        clear templatePad temp;
        
        % apply Gaussian filter to phase correlation data to weight low frequencies
        % more heavily since SNR is higher there
        gaussianFilter = abs(fft2(fspecial('gaussian', [cropHeight cropWidth], p.gaussianFilterSigma)));
        
    end
    
    %% Identify frames to analyze

%     if ~isequal(logPath,0)
%         if length(logFile) == length(dataFile)
%             sifLogData =  importdata([logPath logFile{stack}]);
%             sifLogData = sifLogData(1:numFrames,:);
%         else
%             sifLogData =  importdata([logPath logFile{1}]);
%             sifLogData = sifLogData(frameNum:frameNum+numFrames-1,:);
%             frameNum = frameNum + numFrames;
%         end
%         if channel == '0'
%             selectedFrames = frames;
%         elseif channel == 'G'
%             selectedFrames = find(sifLogData(:,2) == 1);
%         elseif channel == 'R'
%             selectedFrames = find(sifLogData(:,3) == 1);
%         end
%     else
        selectedFrames = frames;
%     end
    
    %% do template matching
    
    hMatchFig = figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
    totalPSFfits = zeros(10000, 6+15+3);
    numPSFfits = 0;
    startTime = tic;
    
    indimall=1;
    
    for a=length(frames):-1:1
        
        if ~logical(sum(frames(a)==selectedFrames))
            continue
        end
        
        data = double(imread([dataPath dataFile{stack}],frames(a),'Info',fileInfo))-darkAvg;
        data = data(ROI(2):ROI(2)+ROI(4)-1, ROI(1):ROI(1)+ROI(3)-1);

 
        % subtract the background and continue
        [bkgndImg,~] = f_waveletBackground(data);
        data = data - bkgndImg;
        
        dataFT = fft2(data,cropHeight,cropWidth);
        maxPeakImg = zeros(cropHeight,cropWidth);
        % matrix PSFLocs stores information about double helices that were
        % found via template matching
        % rows are different matches
        % [xLocation yLocation matchingTemplateNumber matchConfidence];
        PSFLocs = zeros(100,4);
        numPSFLocs = 0;
        for b=1:numTemplates
            % try no prefiltering
            %H = 1;
            % try phase correlation
            %H = 1./(abs(dataFT).*abs(squeeze(templateFT(b,:,:))));
            % try weighted phase correlation (emphasizing low frequency
            % components
            H = gaussianFilter./(abs(dataFT).*abs(squeeze(templateFT(b,:,:))));
            
            % normalize H so it doesn't add any energy to template match
            %H = H / sqrt(sum(abs(H(:)).^2));
            
            peakImg = ifftshift(ifft2(dataFT.*squeeze(templateFT(b,:,:)).*H));
            
            % normalize response of peakImg by dividing by number of pixels in
            % data
            %peakImg = peakImg / (cropHeight*cropWidth);
            maxPeakImg = max(maxPeakImg, peakImg);
            
            % only remember matches that are 3 standard deviations above the
            % mean
            peakThreshold = mean(peakImg(:))+3*std(peakImg(:));
            peakImg(peakImg < peakThreshold) = peakThreshold;
            temp = find(imregionalmax(peakImg));
            % make sure threshold didn't eliminate all peaks and create
            % lots of matches
            if length(temp) < cropHeight*cropWidth/2;
                [tempY, tempX] = ind2sub([cropHeight cropWidth],temp);
                PSFLocs(numPSFLocs+(1:length(temp)),:) = ...
                    [tempX tempY b*ones(length(temp),1) peakImg(temp)];
                numPSFLocs = numPSFLocs+length(temp);
            end
        end
        clear H dataFT peakImg tempX tempY temp
        
        %% filter out extraneous matches due to very strong signals
        
        if numPSFLocs > 0
            % sort location matrix in decending order of confidence
            temp = sortrows(PSFLocs(1:numPSFLocs,:),-4);
            % copy most confident match to list of locations
            PSFLocs(1,:) = temp(1,:);
            numPSFLocs = 1;
            for b=2:size(temp,1)
                % make sure that this candidate location is a minimum distance away
                % from all other candidate locations
                if sum((temp(b,1)-PSFLocs(1:numPSFLocs,1)).^2 + (temp(b,2)-PSFLocs(1:numPSFLocs,2)).^2 >= p.minDistBetweenSMs^2) == numPSFLocs
                    % add it to list of locations
                    numPSFLocs = numPSFLocs + 1;
                    PSFLocs(numPSFLocs,:) = temp(b,:);
                end
            end
        end
        
        totalPSFfits(numPSFfits+1:numPSFfits+numPSFLocs,1:6) = ...
            [frames(a)*ones(numPSFLocs,1) (1:numPSFLocs)' PSFLocs(1:numPSFLocs,:)];
        
        %% output an example image for each threshold level
        %  so that user can pick appropriate threshold later
        
        for b=1:numPSFLocs
            moleThreshold = round(PSFLocs(b,4)*10000);
            moleFileName = ['template ' num2str(PSFLocs(b,3)) ' threshold ' num2str(moleThreshold,'%g') '.png']; %%f6.4 without scaling
            if isempty(dir([outputFilePrefix{stack} moleFileName]))
                % create indices to isolate image of candidate molecule
                [xIdx, yIdx] = meshgrid(PSFLocs(b,1)-p.boxRadius:PSFLocs(b,1)+p.boxRadius, ...
                    PSFLocs(b,2)-p.boxRadius:PSFLocs(b,2)+p.boxRadius);
                % make sure indices are inside ROI
                if min(xIdx(:)) < 1
                    xIdx = xIdx + (1-min(xIdx(:)));
                end
                if max(xIdx(:)) > cropWidth
                    xIdx = xIdx - (max(xIdx(:))-cropWidth);
                end
                if min(yIdx(:)) < 1
                    yIdx = yIdx + (1-min(yIdx(:)));
                end
                if max(yIdx(:)) > cropHeight
                    yIdx = yIdx - (max(yIdx(:))-cropHeight);
                end
                
                % output a picture of the
                img = data(yIdx(:,1),xIdx(1,:));
                img = 1+round(255*(img-min(img(:)))/(max(img(:))-min(img(:))));
%                 imwrite(imresize(ind2rgb(img,hot(256)),3,'nearest'),[outputFilePrefix{stack} moleFileName]);
                imall(indimall).image=img;
                imall(indimall).threshold=moleThreshold;
                imall(indimall).threshold=moleThreshold;
                imall(indimall).template=PSFLocs(b,3);
                indimall=indimall+1;
            end
        end
        numPSFfits = numPSFfits+numPSFLocs;
        
        %%  plot results of template matching and fitting
        if ploton
            set(0,'CurrentFigure',hMatchFig);
            subplot('Position',[0.025 0.025 .9/2 .95],'parent',hMatchFig);
            imagesc(maxPeakImg,[0 3*peakThreshold]);axis image;
            title({'Peaks correspond to likely template matches' ...
                [num2str(numPSFLocs) ' matches found']});

            subplot('Position',[0.525 0.025 .9/2 .95],'parent',hMatchFig);
            imagesc(data);axis image;colormap hot;
            hold on;
            for b=1:numPSFLocs
                plot(PSFLocs(b,1), PSFLocs(b,2), 'o', ...
                    'MarkerSize', 15*PSFLocs(b,4)/peakThreshold, ...
                    'MarkerEdgeColor', templateColors(PSFLocs(b,3),:));
            end
            hold off;
            title({['Frame ' num2str(frames(a)) ': raw data - darkAvg counts'] ...
                ['ROI [xmin ymin width height] = ' mat2str(ROI)]});

            drawnow;
        end
        if indimall>p.maxfits
            break
        end
    end
    elapsedTime = toc(startTime);
    totalPSFfits = totalPSFfits(1:numPSFfits,:);
    clear data bkgnd residual fileInfo maxPeakImg reconstructImg xIdx yIdx temp;
    close(hMatchFig) % closes fitting figure to prevent messiness
    
    fps = length(frames)/elapsedTime
    moleculesPerSec = numPSFfits/elapsedTime
    
    %% output data to external files
    
    textHeader = {'frame number' 'molecule number' 'template x location in ROI (px)' 'template y location in ROI (px)' ...
        'matching template number' 'match confidence (au)' ...
        'amp 1 (counts)' 'amp 2 (counts)' 'x location 1 (px)' 'y location 1 (px)' ...
        'x location 2 (px)' 'y location 2 (px)' 'sigma 1 (px)' 'sigma 2 (px)' 'background mean (counts)' ...
        'total fit error (counts)' 'good fit flag' 'x center (nm)' 'y center (nm)' ...
        'angle (deg)' 'number of photons' 'aberration corrected x location (nm)' ...
        'aberration corrected y location (nm)' 'z location (nm)'};
    % save fit info to MATLAB mat file
    save([outputFilePrefix{stack} 'threshold output.mat']);

    
end

end