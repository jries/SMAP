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

% function [outputFilePrefix] = ...
%     f_fitSMs(dataFile,dataPath,calFile,p.calBeadIdx,templateFile,p.templateFrames,p.peakThreshold,...
%     darkFile,logFile,logPath,p.boxRadius,channel, p.sigmaBounds,p.gaussianFilterSigma,p.minDistBetweenSMs,...
%     p.lobeDistBounds,conversionGain,p.nmPerPixel,EMGain,p.templateLocs,p.ROI)
function [totalPSFfits] = ...
    f_fitSMs(obj,p)

% f_fitSMs is a module in easy_dhpsf that finds the positions of likely 
% DH-PSF profiles by matching to a series of templates generated in 
% f_calDHPSF and prepared in f_calSMidentification. These are then more 
% precisely localized using a double gaussian fit and corrected for drift 
% using the results from f_trackFiducials.

ploton=false;

[dataPath,fi,ex]=fileparts(p.dataFile);
dataPath=[dataPath filesep];
dataFile={[fi ex]};

conversionFactor = p.conversion/p.EMgain;

ampRatioLimit = 0.5;
sigmaRatioLimit = 0.4;

% Options for lsqnonlin
options = optimset('FunValCheck','on','Diagnostics','off','Jacobian','on', 'Display', 'off');
%    'FinDiffType','central','DerivativeCheck','on');



channel='0';
%% ask user for relevant datafiles

outputFilePrefix = cell(1,length(dataFile));
templateFile=p.templateFile;
for stack = 1:length(dataFile)
    
    fileInfo = imfinfo([dataPath dataFile{stack}]);
    numFrames = length(fileInfo);
    frames = 1:numFrames;
    imgHeight = fileInfo(1).Height;
    imgWidth = fileInfo(1).Width;
    
    frames=frames(frames>=p.framerange(1)&frames<=p.framerange(2));
    
    if stack == 1
        
        if strcmp(templateFile(length(templateFile)-2:length(templateFile)),'tif')
            templateInfo = imfinfo(templateFile);
            if templateInfo(1).Height ~= templateInfo(1).Width
                error('Template is not square');
            end
            templateSize = templateInfo(1).Height;
        else
            load(templateFile);
            templateSize = size(template,2);
        end
        
%         temp = inputdlg({'Input frames to process for finding DH-PSFs (ex. [100:199 205 380])'},...
%                          'Input fitting parameters',...
%                          1,...
%                          {['[1:' num2str(length(fileInfo)) ']']});        
% 
%         frames = str2num(temp{1}); %% This only is called for the first iteration of stack =1

    end
    
    %% create output log filenames
    if channel == '0'
        outputFilePrefix{stack} = [dataPath dataFile{stack}(1:length(dataFile{stack})-4) filesep 'molecule fits  ' ...
            datestr(now,'yyyymmdd HHMM') filesep];
        mkdir(outputFilePrefix{stack});
    elseif channel == 'G'
        outputFilePrefix{stack} = [dataPath dataFile{stack}(1:length(dataFile{stack})-4) filesep 'reflected' filesep 'molecule fits  ' ...
            datestr(now,'yyyymmdd HHMM') filesep];
        mkdir(outputFilePrefix{stack});
    elseif channel == 'R'
        outputFilePrefix{stack} = [dataPath dataFile{stack}(1:length(dataFile{stack})-4) filesep 'transmitted' filesep 'molecule fits ' ...
            datestr(now,'yyyymmdd HHMM') filesep];
        mkdir(outputFilePrefix{stack});
    end
    
    if stack == 1
        % Compute darkAvg counts
        
%         if ~isequal(darkFile,0)
%             % Computes average of darkAvg frames for background subtraction
%             darkFileInfo = imfinfo(darkFile);
%             numDarkFrames = length(darkFileInfo);
%             darkAvg = zeros(darkFileInfo(1).Height,darkFileInfo(1).Width);
%             for frame = 1:numDarkFrames
%                 darkAvg = darkAvg + double(imread(darkFile,frame,'Info',darkFileInfo));
%             end
%             darkAvg = darkAvg/numDarkFrames;
%             if ~isequal(size(darkAvg),[imgHeight imgWidth])
%                 warning('Dark count image and data image stack are not the same size. Resizing dark count image...');
%                 darkAvg = imresize(darkAvg,[imgHeight imgWidth]);
%             end
%         else
            darkAvg = p.offsetADU;
%         end
%         clear darkFileInfo;
        
        % Compute average image
        
%         avgImg = zeros(imgHeight,imgWidth);
%         for a = 1:200
%             avgImg = avgImg + double(imread([dataPath dataFile{stack}],a,'Info',fileInfo)) - darkAvg;
%             %    avgImg = avgImg + double(imread([dataPath dataFile],a)) - darkAvg;
%             %Deleted 'Info' -AC 6/21
%         end
%         avgImg = avgImg/200;
        
        % define variables related to the templates
 
        numTemplates = length(p.templateFrames); 
        templateColors = jet(numTemplates);
        
        % make sure p.ROI is an even number of pixels: should also be done in
        % f_calSMidentification
        if mod(p.ROI(3),2)==1
            p.ROI(3) = p.ROI(3)-1;
        end
        if mod(p.ROI(4),2)==1
            p.ROI(4) = p.ROI(4)-1;
        end
        cropWidth = p.ROI(3);
        cropHeight = p.ROI(4);

        %% prepare template for template matching
        
        % pad template to same size as input
        templatePad = zeros(numTemplates,cropHeight,cropWidth);
        templateFT = zeros(numTemplates,cropHeight,cropWidth);
        for a=1:numTemplates
            
            if strcmp(templateFile(length(templateFile)-2:length(templateFile)),'tif')
                templatePad(a,:,:) = padarray(squeeze(template(a,:,:)),...
                    [(cropHeight-size(template,2))/2 ...
                    (cropWidth-size(template,3))/2],min(min(template(a,:,:))));
            else
                templatePad(a,:,:) = padarray(squeeze(template(p.templateFrames(a),:,:)),...
                    [(cropHeight-size(template,2))/2 ...
                    (cropWidth-size(template,3))/2],min(min(template(p.templateFrames(a),:,:))));
            end
            
            % multiplying by conjugate of template in FT domain is equivalent
            % to flipping the template in the real domain
            templateFT(a,:,:) = conj(fft2(squeeze(templatePad(a,:,:))));
        end
        clear templatePad;
        
        % apply Gaussian filter to phase correlation data to weight low frequencies
        % more heavily since SNR is higher there
        gaussianFilter = abs(fft2(fspecial('gaussian', [cropHeight cropWidth], p.gaussianFilterSigma)));
        
        numFrames = length(frames);
        
    end
    
    timecount=tic;
    %% Identify frames to analyze
%     frameNum = 1;
%     
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
    if ploton
    scrsz = get(0,'ScreenSize');
    hSMFits=figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
    end
    totalPSFfits = zeros(20000, 6+18+3);
    numPSFfits = 0;
    startTime = tic;
    logFlag = 0;
    cropWidth = p.ROI(3);
    cropHeight = p.ROI(4);
    bkgndImg = zeros(length(p.ROI(2):p.ROI(2)+p.ROI(4)-1),...
        length(p.ROI(1):p.ROI(1)+p.ROI(3)-1));
    numbkgndImg = 0;
    
    for a=1:numFrames

        if ~logical(sum(frames(a)==selectedFrames))
            continue
        end
        
        data = double(imread([dataPath dataFile{stack}],frames(a),'Info',fileInfo))-darkAvg;
        data = data(p.ROI(2):p.ROI(2)+p.ROI(4)-1, p.ROI(1):p.ROI(1)+p.ROI(3)-1);

        % subtract the background and continue
        bkgndImg_curr = f_waveletBackground(data);
        bkgndImg = bkgndImg + bkgndImg_curr;
        numbkgndImg = numbkgndImg +1;
        data = data - bkgndImg_curr;

        dataFT = fft2(data,cropHeight,cropWidth);
        maxPeakImg = zeros(cropHeight,cropWidth);
        % matrix PSFLocs stores information about double helices that were
        % found via template matching
        % rows are different matches
        % [xLocation yLocation matchingTemplateNumber matchConfidence];
        PSFLocs = zeros(100,4);
        numPSFLocs = 0;
        if length(p.peakThreshold)<numTemplates
            p.peakThreshold(end+1:numTemplates)=p.peakThreshold(1);
        end
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
            
            %threshold = mean(peakImg(:))+p.peakThreshold*std(peakImg(:));
%             peakImg(peakImg < p.peakThreshold(stack,b)) = p.peakThreshold(stack,b);
            peakImg(peakImg < p.peakThreshold(b)) = p.peakThreshold(b);
            
            if isreal(peakImg) && sum(sum(isnan(peakImg)))==0
                temp = find(imregionalmax(peakImg));
            else
                % Write log file
                fileID = fopen('peakImg log.txt','a');
                fprintf(fileID,[datestr(now) '\r\n' dataFile{stack}]);
                fprintf(fileID,'\r\nROI: [%d %d %d %d]',p.ROI);
                fprintf(fileID,'\r\nFrame: %d\r\npeakImg matrix:\r\n',a);
                dlmwrite('peakImg log.txt',peakImg,'-append','delimiter','\t','newline','pc')
                fprintf(fileID,'\r\n********************NEXT********************\r\n\r\n');
                fclose('all');
                
                logFlag = logFlag + 1;
                peakImg(isnan(peakImg)) = 0; %inserted to deal with NaNs -AC 6/22
                temp = find(imregionalmax(real(peakImg)));
            end
            
            
            % make sure threshold didn't eliminate all peaks and create
            % lots of matches
            if length(temp) < cropHeight*cropWidth/2;
                [tempY, tempX] = ind2sub([cropHeight cropWidth],temp);
                PSFLocs(numPSFLocs+(1:length(temp)),:) = ...
                    [tempX tempY b*ones(length(temp),1) peakImg(temp)];
                numPSFLocs = numPSFLocs+length(temp);
            end
        end
        clear H dataFT peakImg
     
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
       
        %% do fitting to extract exact locations of DH-PSFs
        
        % [amp1 amp2 xMean1 yMean1 xMean2 yMean2 sigma1 sigma2 bkgndMean
        %  totalFitError goodFit xCenter yCenter angle numPhotons interlobeDistance amplitude ratio sigma ratio]
        PSFfits = zeros(numPSFLocs, 18);
        % create reconstructed DH-PSF image from fitted data
        reconstructImg = zeros(cropHeight, cropWidth);
        
        for b=1:numPSFLocs
            %% prepare parameters for fine fitting of PSFs
            
            % create indices to use for fitting
            [xIdx, yIdx] = meshgrid(PSFLocs(b,1)-p.boxRadius:PSFLocs(b,1)+p.boxRadius, ...
                PSFLocs(b,2)-p.boxRadius:PSFLocs(b,2)+p.boxRadius);
            % make sure indices are inside p.ROI
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
            
            % compute initial parameters from the location of two spots in
            % the templates
            fitParam(3) = PSFLocs(b,1) + p.templateLocs(PSFLocs(b,3),1)-(templateSize/2+0.5);
            fitParam(4) = PSFLocs(b,2) + p.templateLocs(PSFLocs(b,3),2)-(templateSize/2+0.5);
            fitParam(5) = PSFLocs(b,1) + p.templateLocs(PSFLocs(b,3),3)-(templateSize/2+0.5);
            fitParam(6) = PSFLocs(b,2) + p.templateLocs(PSFLocs(b,3),4)-(templateSize/2+0.5);
            % make sure initial guess lies within the p.ROI: if not, move on
            if fitParam(3)<min(xIdx(:)) || fitParam(5)<min(xIdx(:)) ...
                    || fitParam(3)>max(xIdx(:)) || fitParam(5)>max(xIdx(:)) ...
                    || fitParam(4)<min(yIdx(:)) || fitParam(6)<min(yIdx(:)) ...
                    || fitParam(4)>max(yIdx(:)) || fitParam(6)>max(yIdx(:))
                PSFfits(b,11) = -1000;
                continue;
            end
            %  fitParam(1) = data(round(fitParam(4)),round(fitParam(3)))-bkgndMean;
            %  fitParam(2) = data(round(fitParam(6)),round(fitParam(5)))-bkgndMean;
            fitParam(1) = data(round(fitParam(4)),round(fitParam(3)));
            fitParam(2) = data(round(fitParam(6)),round(fitParam(5)));
            fitParam(7) = mean(p.sigmaBounds);
            fitParam(8) = mean(p.sigmaBounds);
            lowerBound = [0 0 min(xIdx(:)) min(yIdx(:)) min(xIdx(:)) min(yIdx(:)) ...
                p.sigmaBounds(1) p.sigmaBounds(1)];
            %             upperBound = [max(max(data(yIdx(:,1),xIdx(1,:))))-bkgndMean ...
            %                 max(max(data(yIdx(:,1),xIdx(1,:))))-bkgndMean ...
            %                 max(xIdx(:)) max(yIdx(:)) max(xIdx(:)) max(yIdx(:)) ...
            %                 p.sigmaBounds(2) p.sigmaBounds(2)];
            upperBound = [max(max(data(yIdx(:,1),xIdx(1,:)))) ...
                max(max(data(yIdx(:,1),xIdx(1,:)))) ...
                max(xIdx(:)) max(yIdx(:)) max(xIdx(:)) max(yIdx(:)) ...
                p.sigmaBounds(2) p.sigmaBounds(2)];
            
            %% Fit with lsqnonlin
            
            %             [fitParam,temp,residual,exitflag] = lsqnonlin(@(x) ...
            %                 doubleGaussianVector(x,data(yIdx(:,1),xIdx(1,:)),bkgndMean,xIdx,yIdx),...
            %                 fitParam,lowerBound,upperBound,options);
            [fitParam,~,residual,exitflag] = lsqnonlin(@(x) ...
                f_doubleGaussianVector(x,data(yIdx(:,1),xIdx(1,:)),0,xIdx,yIdx),...
                fitParam,lowerBound,upperBound,options);
            
            fittedBkgndMean = mean2(bkgndImg_curr(yIdx(:,1),xIdx(1,:)))*conversionFactor;
            
            PSFfits(b,1:11) = [fitParam fittedBkgndMean sum(abs(residual)) exitflag];
            
            %% compute derived paramters from fine fitting output
            
            % Calculate midpoint between two Gaussian spots
            % shift coordinates relative to entire dataset (not just p.ROI)
            PSFfits(b,3) = PSFfits(b,3) + p.ROI(1)-1;
            PSFfits(b,4) = PSFfits(b,4) + p.ROI(2)-1;
            PSFfits(b,5) = PSFfits(b,5) + p.ROI(1)-1;
            PSFfits(b,6) = PSFfits(b,6) + p.ROI(2)-1;
            % shift coordinates relative to entire dataset (not just p.ROI) and
            % convert from pixels to nm
            PSFfits(b,12) = ((fitParam(3)+fitParam(5))/2 + p.ROI(1)-1)*p.nmPerPixel;
            PSFfits(b,13) = ((fitParam(4)+fitParam(6))/2 + p.ROI(2)-1)*p.nmPerPixel;
            
            % Below is the calculation of the angle of the two lobes.
            % Remember that two vertical lobes is focal plane because camera
            % outputs data that is rotated. Therefore, we want y2>y1 for all
            % angle calculations (so that -90<=angle<=90, and we use swap
            % the use of x and y for the atan2 calculation.
            x1 = fitParam(3);
            x2 = fitParam(5);
            y1 = fitParam(4);
            y2 = fitParam(6);
            % swap if y1>y2
            if (y1 > y2)
                tx = x1; ty = y1;
                x1 = x2; y1 = y2;
                x2 = tx; y2 = ty;
                clear tx ty;
            end
            %Finds the angle
            PSFfits(b,14) = atan2(-(x2-x1),y2-y1) * 180/pi;
            clear x1 x2 y1 y2;
            
            %Below is a way to count the photons. It integrates the box and
            %subtracts the boxarea*offset from the fit. It is inherently flawed
            %if there happens to be bright pixels inside of the fitting region.
            totalCounts = sum(sum(data(yIdx(:,1),xIdx(1,:))));
            %  PSFfits(b,15) = (totalCounts-(2*p.boxRadius+1)^2*bkgndMean)*conversionFactor;  % bkgndMean = 0
            PSFfits(b,15) = totalCounts*conversionFactor;  % bkgndMean = 0

            %The interlobe distance
            lobeDist = sqrt((fitParam(3)-fitParam(5)).^2 + ...
                (fitParam(4)-fitParam(6)).^2);
            PSFfits(b,16) = lobeDist;
            
            %Amplitude Ratio
            ampRatio = abs(fitParam(1) - fitParam(2))/sum(fitParam(1:2));
            PSFfits(b,17) = ampRatio;
            
            % Gaussian width Ratio
            sigmaRatio = abs(fitParam(7) - fitParam(8))/sum(fitParam(7:8));
            PSFfits(b,18) = sigmaRatio;
    
            %% Now evaluate the goodness of the fits
            
            % Conditions for fits (play with these):
            % (1) Amplitude of both lobes > 0
            % (2) All locations x1,y1, x2,y2 lie inside area of small box
            % (3) All sigmas need to be > sigmaBound(1) and < sigmaBound(2)
            % (4) Distance between lobes needs to be > lobeDist(1) pixels and < lobeDist(2) pixels
            % (5) Make sure amplitudes are within 100% of one another
            % (6) Make sure totalFitError/(total number of photons) < 1.05
            
            if exitflag > 0
                % absolute amplitude > 0?
                if fitParam(1)<0 || fitParam(2)<0
                    PSFfits(b,11) = -1001;
                end
                % peaks inside box?
                if fitParam(3)<min(xIdx(:)) || fitParam(5)<min(xIdx(:)) ...
                        || fitParam(3)>max(xIdx(:)) || fitParam(5)>max(xIdx(:)) ...
                        || fitParam(4)<min(yIdx(:)) || fitParam(6)<min(yIdx(:)) ...
                        || fitParam(4)>max(yIdx(:)) || fitParam(6)>max(yIdx(:))
                    PSFfits(b,11) = -1002;
                end
                % absolute sigma size for either lobe within bounds?
                if fitParam(7)<=p.sigmaBounds(1) || fitParam(8)<=p.sigmaBounds(1) ...
                        || fitParam(7)>=p.sigmaBounds(2) || fitParam(8)>=p.sigmaBounds(2)
                    PSFfits(b,11) = -1003;
                end
                % sigma ratio of lobes less than limit?
                if sigmaRatio > sigmaRatioLimit;
                    PSFfits(b,11) = -1004;
                end
                % interlobe distance within bounds?
                if lobeDist < p.lobeDistBounds(1) || lobeDist > p.lobeDistBounds(2)
                    PSFfits(b,11) = -1005;
                end
                % amplitude ratio of lobes less than limit?
                if ampRatio > ampRatioLimit;
                    PSFfits(b,11) = -1006;
                end
                % normalized error within limit?
                if PSFfits(b,10)*conversionFactor/PSFfits(b,15) > 3.0  || ...
                        PSFfits(b,10)*conversionFactor/PSFfits(b,15) < 0.0
                    PSFfits(b,11) = -1007;
                end
                
            end
            
            % if the fit is good, add it to the reconstructed image
            if PSFfits(b,11) > 0
                [xIdx, yIdx] = meshgrid(1:cropWidth,1:cropHeight);
                reconstructImg = reconstructImg + ...
                    fitParam(1).*exp( -((xIdx-fitParam(3)).^2+(yIdx-fitParam(4)).^2.) / (2.*fitParam(7).^2)) ...
                    +fitParam(2).*exp( -((xIdx-fitParam(5)).^2+(yIdx-fitParam(6)).^2.) / (2.*fitParam(8).^2));
            end
        end

        totalPSFfits(numPSFfits+1:numPSFfits+numPSFLocs,6+1:6+18) = PSFfits;
        numPSFfits = numPSFfits+numPSFLocs;
        
        %%  plot results of template matching and fitting
        if ploton
            set(0,'CurrentFigure',hSMFits);
            subplot('Position',[0.025 0.025 .85/3 .95]);
            imagesc(maxPeakImg,[0 3*min(p.peakThreshold(stack,:))]);axis image;
            title({'Peaks correspond to likely template matches' ...
                [num2str(numPSFLocs) ' matches found']});

            subplot('Position',[0.075+.85/3 0.025 .85/3 .95]);
            imagesc(data);axis image;colormap hot;
            hold on;
            for b=1:numPSFLocs
                %plot(PSFLocs(b,1), PSFLocs(b,2), 'o', ...
                %    'MarkerSize', 15*PSFLocs(b,4)/p.peakThreshold(b), ...
                %    'MarkerEdgeColor', templateColors(PSFLocs(b,3),:));
                plot(PSFLocs(b,1), PSFLocs(b,2), 'o', ...
                    'MarkerSize', 15*PSFLocs(b,4)/p.peakThreshold(PSFLocs(b,3)), ...
                    'MarkerEdgeColor', templateColors(PSFLocs(b,3),:));  
                %Jonas: stack removed
                %Changed p.peakThreshold(b) to p.peakThreshold(PSFLocs(b,3)) because b
                %does not seem to logically correspond to the correct template, whereas the third
                %column of PSFLocs is defined as corresponding to a specific
                %template. Furthermore, whenever numPSFLocs > length(p.peakThreshold)
                %there is an error. -AC 6/22
            end
            hold off;
            title({['Frame ' num2str(frames(a)) ': raw data - darkAvg counts'] ...
                ['p.ROI [xmin ymin width height] = ' mat2str(p.ROI)]});

            subplot('Position',[0.125+2*.85/3 0.025 .85/3 .95]);
    %         imagesc(reconstructImg+bkgndMean,[min(data(:)) max(data(:))]);axis image;
            imagesc(reconstructImg,[min(data(:)) max(data(:))]);axis image;
            title({'Image reconstructed from fitted matches' ...
                [num2str(sum(PSFfits(:,11)>0)) ' successful fits']});

            drawnow;
        end
        %M(a) = getframe(h);
        %imwrite(frame2im(M),'test_output.tif','tif','WriteMode','append');
        if toc(timecount)>1
        obj.status(['frame processed: ' num2str(a) ' of ' num2str(numFrames)]); drawnow;
        timecount=tic;
        end
    end
    elapsedTime = toc(startTime);
    totalPSFfits = totalPSFfits(1:numPSFfits,:);
    clear data bkgnd residual fileInfo maxPeakImg reconstructImg xIdx yIdx temp;
    %     clear data bkgnd residual fileInfo maxPeakImg reconstructImg templateFT xIdx yIdx temp;
    
    if logFlag ~= 0
        sprintf('There were %d frames in which the peakImg matrix contained complex numbers or NaNs. See log file "peakImg log.txt" for more details',logFlag)
        logFlag = 0;
    end
    if ploton
    close(hSMFits)
    end
    fps = length(frames)/elapsedTime
    moleculesPerSec = numPSFfits/elapsedTime
    %movie2avi(M,'output_v1.avi','fps',10,'Compression','FFDS');
    
    %% Measure the Gaussian Laser Intensity Distribution
%     bkgndImg_avg = bkgndImg/numbkgndImg;
    
%     [laser_x_nm, laser_y_nm ,sigma_x_nm, sigma_y_nm, theta, peakIntensity, waist]...
%         = EstimateGaussianLaserProfile...
%         (bkgndImg_avg, FOWmask, p.nmPerPixel, powerAtObjective, p.ROI);
    
    
    save([outputFilePrefix{stack} 'molecule fits.mat']);
    
    %% Translate angle into corrected x,y positions and z position
    
    load(p.calFile);
    totalPSFfits(:,25) = totalPSFfits(:,18) ...
        - interp1(squeeze(meanAngles(1,p.calBeadIdx,goodFit_forward)),...
        squeeze(meanX(1,p.calBeadIdx,goodFit_forward)),totalPSFfits(:,20),'spline');
    totalPSFfits(:,26) = totalPSFfits(:,19) ...
        - interp1(squeeze(meanAngles(1,p.calBeadIdx,goodFit_forward)),...
        squeeze(meanY(1,p.calBeadIdx,goodFit_forward)),totalPSFfits(:,20),'spline');
    totalPSFfits(:,27) = interp1(squeeze(meanAngles(1,p.calBeadIdx,goodFit_forward)),...
        z(goodFit_forward),totalPSFfits(:,20),'spline');
    
    %% output data to external file
   
    save([outputFilePrefix{stack} 'molecule fits.mat']);

    
end
end
