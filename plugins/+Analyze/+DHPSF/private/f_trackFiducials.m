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

function [outputFilePrefix] = f_trackFiducials(dataFile,dataPath,calFile,calBeadIdx,templateFile,templateFrames,...
                     darkFile,logFile,logPath,boxRadius,channel, gaussianFilterSigma,minDistBetweenSMs,...
                     lobeDistBounds,conversionGain,nmPerPixel,EMGain,templateLocs,sigmaBounds)
% f_trackFiducials is a module in easy_dhpsf that extracts the position of
% one or more fiducial beads in an image stack. This position is then used 
% to correct the fit localizations from that image stack.

 % sigmaBounds = [1.0 3.0] * 125.78 / nmPerPixel;
    
conversionFactor = conversionGain/EMGain;

% sigmaBounds = [1.0 3.0];       % sets [min max] allowed sigma for double Gaussian fit (in units of pixels)
% lobeDistBounds = [6.0 10.0];   % sets [min max] allowed interlobe distance for double Gaussian fit (in units of pixels)
ampRatioLimit = 0.4;            % sets maximum allowed amplitude ratio
sigmaRatioLimit = 0.3;          % sets maximum allowed sigma ratio

scrsz = get(0,'ScreenSize');
numSyncFrames = 25;
% Options for lsqnonlin
options = optimset('FunValCheck','on','Diagnostics','off','Jacobian','on', 'Display', 'off');
%    'FinDiffType','central','DerivativeCheck','on');


outputFilePrefix = cell(1,length(dataFile));

for stack = 1:length(dataFile)
    

    
    dataFileInfo = imfinfo([dataPath dataFile{stack}]);
    numFrames = length(dataFileInfo);
    frames = 1:numFrames;
    imgHeight = dataFileInfo(1).Height;
    imgWidth = dataFileInfo(1).Width;
    
    if stack == 1
  
        load(calFile)
          
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

        peakThreshold = 0.01*ones(length(templateFrames),1);
        
        % Construct a questdlg with three options
%         dlg_title = 'Please Verify';
%         prompt = {  'dataFile = ' [dataFile{stack} ' to ' dataFile{length(dataFile)}], ...
%             'darkFile = ' darkFile, ...
%             'logFile = ' [logFile{stack} ' to ' logFile{length(logFile)}]
%             };
%         def =       { 'Yes'  };
%         questiondialog = questdlg(prompt,dlg_title, def);
%         % Handle response
%         switch questiondialog
%             case 'Yes'
%             case 'No'
%                 error('User cancelled the program');
%             case 'Cancel'
%                 error('User cancelled the program');
%         end
        
    end
    if channel == '0'
        outputFilePrefix{stack} = [dataPath dataFile{stack}(1:length(dataFile{stack})-4) '\fiduciaries ' ...
        datestr(now,'yyyymmdd HHMM') '\'];
        mkdir(outputFilePrefix{stack});
    
    elseif channel == 'G'
        outputFilePrefix{stack} = [dataPath dataFile{stack}(1:length(dataFile{stack})-4) '\reflected\fiduciaries ' ...
            datestr(now,'yyyymmdd HHMM') '\'];
        mkdir(outputFilePrefix{stack});
        
    elseif channel == 'R'
        outputFilePrefix{stack} = [dataPath dataFile{stack}(1:length(dataFile{stack})-4) '\transmitted\fiduciaries ' ...
            datestr(now,'yyyymmdd HHMM') '\'];
        mkdir(outputFilePrefix{stack});
    
    end
    %     outputFilePrefix{stack} = [dataPath dataFile{stack}(1:length(dataFile{stack})-4) '\fiduciaries ' ...
%         datestr(now,'yyyymmdd HHMM') '\'];
%     mkdir(outputFilePrefix{stack});
    
    
    if stack == 1
        
        % Compute dark counts
        
        if ~isequal(darkFile,0)
            % Computes average of dark frames for background subtraction
            darkFileInfo = imfinfo(darkFile);
            numDarkFrames = length(darkFileInfo);
            darkAvg = zeros(darkFileInfo(1).Height,darkFileInfo(1).Width);
            for frame = 1:numDarkFrames
                darkAvg = darkAvg + double(imread(darkFile,frame,'Info',darkFileInfo));
            end
            darkAvg = darkAvg/numDarkFrames;
            if ~isequal(size(darkAvg),[imgHeight imgWidth])
                warning('Dark count image and data image stack are not the same size. Resizing dark count image...');
                darkAvg = imresize(darkAvg,[imgHeight imgWidth]);
            end
        else
            darkAvg = 0;
        end
        clear darkFileInfo;
        
        %% Compute average image
        dataAvg = zeros(imgHeight,imgWidth);
        for frame = 1:200
            % for frame = 1:(numFrames/10)
            dataAvg = dataAvg + double(imread([dataPath dataFile{stack}],frame,'Info',dataFileInfo)) - darkAvg;
        end
        % dataAvg = dataAvg/(numFrames/10);
        dataAvg = dataAvg/200;
        
        % definitions related to templates
        numTemplates = length(templateFrames);
        %% Pick region of interest for analysis
        
        hROI = figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
        imagesc(dataAvg);axis image;colormap hot;
        if channel == 'G' || channel == '0'
            ROI = imrect(gca,[1 1 270 270]);
        elseif channel == 'R'
            ROI = imrect(gca,[243 243 270 270]);
        end
        
        title({'Double-click to choose region of interest for PSF extraction' ...
            mat2str(ROI.getPosition)});
        addNewPositionCallback(ROI,@(p) title({'Double-click to choose region of interest for PSF extraction' ...
            ['[xmin ymin width height] = ' mat2str(p,3)]}));
 %       make sure rectangle stays within image bounds
        fcn = makeConstrainToRectFcn('imrect',get(gca,'XLim'),get(gca,'YLim'));
        setPositionConstraintFcn(ROI,fcn);
        ROI = round(wait(ROI));
        
        % if odd, rounds down to nearest even number for width and height 
        % (must be divisible by two for padding step)
        ROI(3) = 2*floor(ROI(3)/2); 
        ROI(4) = 2*floor(ROI(4)/2);
      
        imgHeight = ROI(4);
        imgWidth = ROI(3);
        cropHeight = imgHeight;
        cropWidth = imgWidth;
        close(hROI)
        %% Ask user for molecule location(s)

        % Plots the avg image .tif
        hMLocs = figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
        imagesc(dataAvg(ROI(2):ROI(2)+ROI(4)-1, ...
            ROI(1):ROI(1)+ROI(3)-1));axis image;colorbar;colormap hot;
        title('Use LMB to select fiducials. Hit enter to stop or use RMB to select the final fiducial.');
        hold on;
        % User will click on molecules in image and then text will be imposed on the
        % image corresponding to the molecule number.
        % pointstofit is a n by 2 array -- x y values for each molecule
        % hold on;
        moleLocs = [];
        n = 0;
        % Loop, picking up the points.
        but = 1;
        while but == 1
            [xi,yi,but] = ginput(1);
            if isempty(xi)
                break
            end
            n = n+1;
            text(xi,yi,num2str(n),'color','white','fontsize',13,'fontweight','bold');
            moleLocs(n,:) = round([xi yi]);
        end
        hold off;
        
        %moleLocs = round(ginput);
        numMoles = size(moleLocs,1);
        
        % moleMask = ones(imgHeight,imgWidth);
        % blackoutbox = 26;
        % for a = 1:numMoles
        %     moleMask(moleLocs(a,2)-blackoutbox:moleLocs(a,2)+blackoutbox,...
        %         moleLocs(a,1)-blackoutbox:moleLocs(a,1)+blackoutbox) = 0;
        % end
        % imagesc(dataAvg.*moleMask);
        % axis square
        
        % save the image with the numbered molecules for future reference
        saveas(hMLocs,[outputFilePrefix{stack} 'molecule map.png']);
        close(hMLocs)
        
        %% user picks background levels
%         figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
%         imagesc(dataAvg(ROI(2):ROI(2)+ROI(4)-1, ...
%             ROI(1):ROI(1)+ROI(3)-1),[0 300]);
%         axis image;colorbar;colormap hot;
%         title('Measure [min max] for Gaussian Laser Background');
%         dlg_title = 'Please Enter Thresholds for Gaussian Laser Background Fitting';
%         prompt = {  'Lower Limit:', ...
%             'Upper Limit:' ...
%             };
%         def = {     '20', ...
%             '100' ...
%             };
%         num_lines = 1;
%         inputdialog = inputdlg(prompt,dlg_title,num_lines,def);
%         threshold = [str2double(inputdialog{1}) str2double(inputdialog{2})];
%         
        
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
                templatePad(a,:,:) = padarray(squeeze(template(templateFrames(a),:,:)),...
                    [(cropHeight-size(template,2))/2 ...
                    (cropWidth-size(template,3))/2],min(min(template(templateFrames(a),:,:))));
            end
            
            % multiplying by conjugate of template in FT domain is equivalent
            % to flipping the template in the real domain
            templateFT(a,:,:) = conj(fft2(squeeze(templatePad(a,:,:))));
        end
        clear templatePad;
        
        % apply Gaussian filter to phase correlation data to weight low frequencies
        % more heavily since SNR is higher there
        gaussianFilter = abs(fft2(fspecial('gaussian', [cropHeight cropWidth], gaussianFilterSigma)));
        
        
    end
    
    %% Fit chosen molecules throughout entire image
    hMFit = figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
    
    % [frameNum moleNum amp1 amp2 xMean1 yMean1 xMean2 yMean2 sigma1 sigma2
    %  bkgndMean totalFitError goodFit xCenter yCenter angle numPhotons ...
    %  xAberrCorrected yAberCorrected zLocation]
    PSFfits = zeros(numFrames*numMoles, 23);
    numPSFfits = 0;
    bkgndFits = zeros(numFrames,8);
    startTime = tic;
    %     meanData = zeros(200,1)
    %     meanBGsubData = zeros(200,1)
    imgHeight = ROI(4);
    imgWidth = ROI(3);
    
    %% Identify frames to analyze
    
    if ~isequal(logPath,0)
        sifLogData =  importdata([logPath logFile{stack}]);
        sifLogData = sifLogData(1:numFrames,:);
        if channel == 'G'
            selectedFrames = find(sifLogData(:,2) == 1);
        elseif channel == 'R'
            selectedFrames = find(sifLogData(:,3) == 1);
        end
    else
        selectedFrames = frames;
    end

    for a=1:numFrames
        
        if ~logical(sum(frames(a)==selectedFrames))
            continue
        end
        
        data = double(imread([dataPath dataFile{stack}],a,'Info',dataFileInfo))-darkAvg;
        data = data(ROI(2):ROI(2)+ROI(4)-1, ...
            ROI(1):ROI(1)+ROI(3)-1);
        
        bkgndImg = f_waveletBackground(data);
        data = data - bkgndImg;
        
%         %% compute background image for each frame
%         
%         %    bkgnd = data(bkgndROI(2):bkgndROI(2)+bkgndROI(4)-1, ...
%         %        bkgndROI(1):bkgndROI(1)+bkgndROI(3)-1);
%         %    bkgndMean = mean(bkgnd(:));
%         [xIdx yIdx] = meshgrid(1:size(data,2), 1:size(data,1));
%         if ~exist('bkgndFit','var')
%             bkgndFit = [median(data(:))-min(data(:)) ...
%                 size(data,2)/2 size(data,1)/2 ...
%                 size(data,2)/4 size(data,1)/4 0 min(data(:))];
%         end
%         
%         % Fit with lsqnonlin
%         lowerBound = [0 size(data,2)/5 size(data,1)/5 ...
%             size(data,2)/5 size(data,1)/5 -360 min(data(:))];
%         upperBound = [max(data(:))-min(data(:)) size(data,2)*3/4 size(data,1)*3/4 ...
%             size(data,2)*4/5 size(data,1)*4/5 360 max(data(:))];
%         bkgndFit = lsqnonlin(@(x) singleGaussianRotOffset(x,data,xIdx,yIdx,threshold),...
%             bkgndFit,lowerBound,upperBound);
%         
%         sigma_x = bkgndFit(4);
%         sigma_y = bkgndFit(5);
%         theta = bkgndFit(6)*pi/180;
%         A = cos(theta)^2/2/sigma_x^2 + sin(theta)^2/2/sigma_y^2;
%         B = -sin(2*theta)/4/sigma_x^2 + sin(2*theta)/4/sigma_y^2 ;
%         C = sin(theta)^2/2/sigma_x^2 + cos(theta)^2/2/sigma_y^2;
%         
%         bkgndImg = bkgndFit(1)*exp( -(A*(xIdx-bkgndFit(2)).^2 + 2*B*(xIdx-bkgndFit(2)).*(yIdx-bkgndFit(3)) + C*(yIdx-bkgndFit(3)).^2)) ...
%             +bkgndFit(7);
%         %imwrite(uint16(bkgndImg),[bkgndFile num2str(a) '.tif'],'tif');
%         bkgndFits(a,:) = [a bkgndFit];
        
        
        %% do fitting to extract exact locations of DH-PSFs
        %         meanData(a) = mean2(data(100:126,50:76))*conversionFactor;
%         bkgndImg = 0;
%         data = data - bkgndImg;
        
        %         meanBGsubData(a) = mean2(data(100:126,50:76))*conversionFactor;
        %     imagesc(data,[0 200])
        %     drawnow
        %     pause(2)
%         bkgndMean = 0;
        
        % create reconstructed DH-PSF image from fitted data
        reconstructImg = zeros(imgHeight, imgWidth);
        
        for b=1:numMoles
            
            rowIdx = (a-1)*numMoles+b;
            
            
            dataFT = fft2(data,cropHeight,cropWidth);
            maxPeakImg = zeros(cropHeight,cropWidth);
            % matrix PSFLocs stores information about double helices that were
            % found via template matching
            % rows are different matches
            % [xLocation yLocation matchingTemplateNumber matchConfidence];
            PSFLocs = zeros(100,4);
            numPSFLocs = 0;
            for c=1:numTemplates
                % try no prefiltering
                %H = 1;
                % try phase correlation
                %H = 1./(abs(dataFT).*abs(squeeze(templateFT(b,:,:))));
                % try weighted phase correlation (emphasizing low frequency
                % components
                H = gaussianFilter./(abs(dataFT).*abs(squeeze(templateFT(c,:,:))));
                
                % normalize H so it doesn't add any energy to template match
                %H = H / sqrt(sum(abs(H(:)).^2));
                
                peakImg = ifftshift(ifft2(dataFT.*squeeze(templateFT(c,:,:)).*H));
                
                % normalize response of peakImg by dividing by number of pixels in
                % data
                %peakImg = peakImg / (cropHeight*cropWidth);
                maxPeakImg = max(maxPeakImg, peakImg);
                
                %threshold = mean(peakImg(:))+peakThreshold*std(peakImg(:));
                peakImg(peakImg < peakThreshold(c)) = peakThreshold(c);
                
                if isreal(peakImg) && sum(sum(isnan(peakImg)))==0
                    temp = find(imregionalmax(peakImg));
                else
                    % Write log file
                    fileID = fopen('peakImg log.txt','a');
                    fprintf(fileID,[datestr(now) '\r\n' dataFile{stack}]);
                    fprintf(fileID,'\r\nROI: [%d %d %d %d]',ROI);
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
                    [tempY tempX] = ind2sub([cropHeight cropWidth],temp);
                    PSFLocs(numPSFLocs+(1:length(temp)),:) = ...
                        [tempX tempY c*ones(length(temp),1) peakImg(temp)];
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
                for c=2:size(temp,1)
                    % make sure that this candidate location is a minimum distance away
                    % from all other candidate locations
                    if sum((temp(c,1)-PSFLocs(1:numPSFLocs,1)).^2 + (temp(c,2)-PSFLocs(1:numPSFLocs,2)).^2 >= minDistBetweenSMs^2) == numPSFLocs
                        % add it to list of locations
                        numPSFLocs = numPSFLocs + 1;
                        PSFLocs(numPSFLocs,:) = temp(c,:);
                    end
                end
            end
            
            %         totalPSFfits(numPSFfits+1:numPSFfits+numPSFLocs,1:6) = ...
            %             [frames(a)*ones(numPSFLocs,1) (1:numPSFLocs)' PSFLocs(1:numPSFLocs,:)];
            
            if numPSFLocs == 1
                % create indices to use for fitting
                [xIdx yIdx] = meshgrid(moleLocs(b,1)-boxRadius:moleLocs(b,1)+boxRadius, ...
                    moleLocs(b,2)-boxRadius:moleLocs(b,2)+boxRadius);
                % make sure indices are inside image
                if min(xIdx(:)) < 1
                    xIdx = xIdx + (1-min(xIdx(:)));
                end
                if max(xIdx(:)) > imgWidth
                    xIdx = xIdx - (max(xIdx(:))-imgWidth);
                end
                if min(yIdx(:)) < 1
                    yIdx = yIdx + (1-min(yIdx(:)));
                end
                if max(yIdx(:)) > imgHeight
                    yIdx = yIdx - (max(yIdx(:))-imgHeight);
                end
                
                %             % if fit failed or if this is the first frame, use two largest
                %             % peaks in the box as initial fitting parameters
                %             % otherwise, use parameters from previous fit
                %             %         if a == 1
                %
                %             % find two largest peaks in box
                %             [tempY tempX] = ind2sub([2*boxRadius+1 2*boxRadius+1], ...
                %                 find(imregionalmax(data(yIdx(:,1),xIdx(1,:)))));
                %             tempX = tempX + min(xIdx(:))-1;
                %             tempY = tempY + min(yIdx(:))-1;
                %             temp = sortrows([tempX tempY data(sub2ind([imgHeight imgWidth], ...
                %                 tempY,tempX))],-3);
                
                
                % set initial fitting parameters
                % [amp1 amp2 xMean1 yMean1 xMean2 yMean2 sigma1 sigma2]
                %             fitParam(3) = temp(1,1);
                %             fitParam(4) = temp(1,2);
                %             fitParam(5) = temp(2,1);
                %             fitParam(6) = temp(2,2);
                %             fitParam(1) = temp(1,3);
                %             fitParam(2) = temp(2,3);
                %             lobeDist = sqrt((temp(1,1)-temp(2,1)).^2 + ...
                %                 (temp(1,2)-temp(2,2)).^2)
                %             if  lobeDist < lobeDistBounds(1) || lobeDist > lobeDistBounds(2)
                %                 fitParam(3) = temp(1,1);
                %                 fitParam(4) = temp(1,2);
                %                 fitParam(5) = temp(3,1);
                %                 fitParam(6) = temp(3,2);
                %                 fitParam(1) = temp(1,3);
                %                 fitParam(2) = temp(3,3);
                %             end
                %         elseif PSFfits(rowIdx-numMoles,13) < 0
                %             i = 2;
                %             while PSFfits(rowIdx-i*numMoles,13) < 0
                %                 i = i+1;
                %             end
                %             fitParam(1:6) = PSFfits(rowIdx-i*numMoles,3:8);
                %
                %         else
                %             % set initial parameters equal to parameters from previous fit
                %             fitParam(1:6) = PSFfits(rowIdx-numMoles,3:8);
                %         end
                
                %             fitParam(7) = 1.8;
                %             fitParam(8) = 1.8;
                %             lowerBound = [0 0 min(xIdx(:)) min(yIdx(:)) min(xIdx(:)) min(yIdx(:)) ...
                %                 sigmaBounds(1) sigmaBounds(1)];
                %             upperBound = [max(max(data(yIdx(:,1),xIdx(1,:))))-bkgndMean ...
                %                 max(max(data(yIdx(:,1),xIdx(1,:))))-bkgndMean ...
                %                 max(xIdx(:)) max(yIdx(:)) max(xIdx(:)) max(yIdx(:)) ...
                %                 sigmaBounds(2) sigmaBounds(2)];
                
                
                
                % compute initial parameters from the location of two spots in
                % the templates
                % [amp1 amp2 xMean1 yMean1 xMean2 yMean2 sigma1 sigma2]
                fitParam(3) = PSFLocs(b,1) + templateLocs(PSFLocs(b,3),1)-(templateSize/2+0.5);
                fitParam(4) = PSFLocs(b,2) + templateLocs(PSFLocs(b,3),2)-(templateSize/2+0.5);
                fitParam(5) = PSFLocs(b,1) + templateLocs(PSFLocs(b,3),3)-(templateSize/2+0.5);
                fitParam(6) = PSFLocs(b,2) + templateLocs(PSFLocs(b,3),4)-(templateSize/2+0.5);
                % make sure initial guess lies within the ROI: if not, move on
                if fitParam(3)<min(xIdx(:)) || fitParam(5)<min(xIdx(:)) ...
                        || fitParam(3)>max(xIdx(:)) || fitParam(5)>max(xIdx(:)) ...
                        || fitParam(4)<min(yIdx(:)) || fitParam(6)<min(yIdx(:)) ...
                        || fitParam(4)>max(yIdx(:)) || fitParam(6)>max(yIdx(:))
                    PSFfits(rowIdx,13) = -1000;
                    PSFfits(rowIdx,1:2) = [a b];
                    continue;
                end
                fitParam(1) = data(round(fitParam(4)),round(fitParam(3)));
                fitParam(2) = data(round(fitParam(6)),round(fitParam(5)));
                fitParam(7) = 1.8;
                fitParam(8) = 1.8;
                lowerBound = [0 0 min(xIdx(:)) min(yIdx(:)) min(xIdx(:)) min(yIdx(:)) ...
                    sigmaBounds(1) sigmaBounds(1)];
                upperBound = [max(max(data(yIdx(:,1),xIdx(1,:)))) ...
                    max(max(data(yIdx(:,1),xIdx(1,:)))) ...
                    max(xIdx(:)) max(yIdx(:)) max(xIdx(:)) max(yIdx(:)) ...
                    sigmaBounds(2) sigmaBounds(2)];
                
                
                
                %% Fit with lsqnonlin
                [fitParam,temp,residual,exitflag] = lsqnonlin(@(x) ...
                    f_doubleGaussianVector(x,data(yIdx(:,1),xIdx(1,:)),0,xIdx,yIdx),...
                    fitParam,lowerBound,upperBound,options);
                PSFfits(rowIdx,1:13) = [a b fitParam 0 sum(abs(residual)) exitflag];
                
                % Calculate midpoint between two Gaussian spots
                % convert from pixels to nm
                PSFfits(rowIdx,14) = ((fitParam(3)+fitParam(5))/2)*nmPerPixel;
                PSFfits(rowIdx,15) = ((fitParam(4)+fitParam(6))/2)*nmPerPixel;
                
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
                PSFfits(rowIdx,16) = atan2(-(x2-x1),y2-y1) * 180/pi;
                clear x1 x2 y1 y2;
                
                %Below is a way to count the photons. It integrates the box and
                %subtracts the boxarea*offset from the fit. It is inherently flawed
                %if there happens to be bright pixels inside of the fitting region.
                totalCounts = sum(sum(data(yIdx(:,1),xIdx(1,:))));
                PSFfits(rowIdx,17) = totalCounts*conversionFactor;
                
                %The interlobe distance
                lobeDist = sqrt((fitParam(3)-fitParam(5)).^2 + ...
                    (fitParam(4)-fitParam(6)).^2);
                PSFfits(rowIdx,18) = lobeDist;
                
                %Amplitude Ratio
                ampRatio = abs(fitParam(1) - fitParam(2))/sum(fitParam(1:2));
                PSFfits(b,19) = ampRatio;
                
                % Gaussian width Ratio
                simgaRatio = abs(fitParam(7) - fitParam(8))/sum(fitParam(7:8));
                PSFfits(b,20) = simgaRatio;
                
                
                %% Now evaluate the fits
                % Conditions for fits (play with these):
                % (1) Amplitude of both lobes > 0
                % (2) All locations x1,y1, x2,y2 lie inside area of small box
                % (3) All sigmas need to be > sigmaBound(1) and < sigmaBound(2)
                % (4) Distance between lobes needs to be > lobeDist(1) pixels and < lobeDist(2) pixels
                % (5) Make sure amplitudes are within 100% of one another
                % (6) Make sure totalFitError/(total number of photons) < 1.05
                
                if exitflag > 0
                    if fitParam(1)<0 || fitParam(2)<0
                        PSFfits(rowIdx,13) = -1001;
                    end
                    if fitParam(3)<min(xIdx(:)) || fitParam(5)<min(xIdx(:)) ...
                            || fitParam(3)>max(xIdx(:)) || fitParam(5)>max(xIdx(:)) ...
                            || fitParam(4)<min(yIdx(:)) || fitParam(6)<min(yIdx(:)) ...
                            || fitParam(4)>max(yIdx(:)) || fitParam(6)>max(yIdx(:))
                        PSFfits(rowIdx,13) = -1002;
                    end
                    if fitParam(7)<=sigmaBounds(1) || fitParam(8)<=sigmaBounds(1) ...
                            || fitParam(7)>=sigmaBounds(2) || fitParam(8)>=sigmaBounds(2)
                        PSFfits(rowIdx,13) = -1003;
                    end
                    if simgaRatio > sigmaRatioLimit;
                        PSFfits(rowIdx,13) = -1004;
                    end
                    if lobeDist < lobeDistBounds(1) || lobeDist > lobeDistBounds(2)
                        PSFfits(rowIdx,13) = -1005;
                    end
                    if ampRatio > ampRatioLimit;
                        PSFfits(rowIdx,13) = -1006;
                    end
%                     if PSFfits(rowIdx,12)*conversionFactor/PSFfits(rowIdx,12) > 3.0  || ...
%                             PSFfits(rowIdx,12)*conversionFactor/PSFfits(rowIdx,12) < 0.0
%                         PSFfits(rowIdx,13) = -1007;
%                     end
                    
                end
                
                % if fit was successful, use the computed center location as center
                % of box for next iteration
                if PSFfits(rowIdx,13) > 0
                    moleLocs(b,:) = round(PSFfits(rowIdx,14:15)/nmPerPixel);
                end
                
                % plot image reconstruction so that fits can be checked
                [xIdx yIdx] = meshgrid(1:imgWidth,1:imgHeight);
                reconstructImg = reconstructImg + ...
                    fitParam(1).*exp( -((xIdx-fitParam(3)).^2+(yIdx-fitParam(4)).^2.) / (2.*fitParam(7).^2)) ...
                    +fitParam(2).*exp( -((xIdx-fitParam(5)).^2+(yIdx-fitParam(6)).^2.) / (2.*fitParam(8).^2));
            end
            PSFfits(rowIdx,1:2) = [a b];
        end
        
        %%  plot results of template matching and fitting
        set(0,'CurrentFigure',hMLocs);
        subplot('Position',[0.025 0.025 .9/2 .95]);
        imagesc(data);axis image;colormap hot;
        title(['Frame ' num2str(a) ': raw data - dark counts']);
        
        subplot('Position',[0.525 0.025 .9/2 .95]);
        imagesc(reconstructImg,[min(data(:)) max(data(:))]);axis image;
        title('Image reconstructed from fitted matches');
        
        drawnow;

    end
    elapsedTime = toc(startTime);
    clear data dataFileInfo residual dataAvg reconstructImg xIdx yIdx temp;
    fps = numFrames/elapsedTime
    moleculesPerSec = numFrames*numMoles/elapsedTime
    close(hMLocs)
    %% Translate angle into corrected x,y positions and z position
    
    load(calFile);%squeeze(meanAngles(1,calBeadIdx,goodFit_forward)
    PSFfits(:,21) = PSFfits(:,14) - interp1(squeeze(meanAngles(1,calBeadIdx,goodFit_forward)),squeeze(meanX(1,calBeadIdx,goodFit_forward)),PSFfits(:,16),'spline');
    PSFfits(:,22) = PSFfits(:,15) - interp1(squeeze(meanAngles(1,calBeadIdx,goodFit_forward)),squeeze(meanY(1,calBeadIdx,goodFit_forward)),PSFfits(:,16),'spline');
    PSFfits(:,23) = interp1(squeeze(meanAngles(1,calBeadIdx,goodFit_forward)),z(goodFit_forward),PSFfits(:,16),'spline');
    
    %% Output raw fit data
    
%     textHeader = {'frame number' 'molecule number' ...
%         'amp 1 (counts)' 'amp 2 (counts)' 'x mean 1 (px)' 'y mean 1 (px)' ...
%         'x mean 2 (px)' 'y mean 2 (px)' 'sigma 1 (px)' 'sigma 2 (px)' 'background mean (counts)' ...
%         'total fit error (counts)' 'good fit flag' 'x center (nm)' 'y center (nm)' ...
%         'angle (deg)' 'number of photons' 'interlobe distance' 'amplitude ratio' 'sigma ratio' 'aberration corrected x location (nm)' ...
%         'aberration corrected y location (nm)' 'z location (nm)'};
    % save fit info to MATLAB mat file
    save([outputFilePrefix{stack} 'raw fits.mat']);
%     % write fitted data to Excel spreadsheet
%     xlswrite([outputFilePrefix{stack} 'raw fits.xlsx'], [textHeader; num2cell(PSFfits)], 'PSF fits');
%     xlswrite([outputFilePrefix{stack} 'raw fits.xlsx'], {'Data filename:' [dataPath dataFile{stack}]; ...
%         'Dark count filename:' darkFile; ...
%         'Background region of interest' mat2str(ROI)}, 'fitting settings');
    
    
    %% compute movement of fiduciaries
    
    devX = zeros(numFrames,numMoles);
    devY = zeros(numFrames,numMoles);
    devZ = zeros(numFrames,numMoles);
    goodFitFlag = zeros(numFrames,numMoles);
    numPhotons = zeros(numFrames,numMoles);
    avgDevX = zeros(numFrames,1);
    avgDevY = zeros(numFrames,1);
    avgDevZ = zeros(numFrames,1);
    numValidFits = zeros(numFrames,1);
    
%     textHeader = {'frame number' 'deviation in x (nm)' 'deviation in y (nm)' ...
%         'deviation in z (nm)' 'good fit flag' 'number of photons'};
    
    syncFrames = zeros(1,numSyncFrames);
    lastGoodFrame = numFrames;
    for a = 1:numSyncFrames
        while sum(PSFfits(PSFfits(PSFfits(:,1)==lastGoodFrame,13)>0,2)) ~= numMoles/2*(1+numMoles)
            lastGoodFrame = lastGoodFrame - 1;
            if lastGoodFrame < 0
                error('No good frames: check dataset and fit flags')
            end
        end
        syncFrames(a)=lastGoodFrame;
        lastGoodFrame = lastGoodFrame - 1;
    end
    
    for molecule = 1:numMoles
        %% extract fitting parameters for this molecule
        moleculeFitParam = PSFfits(PSFfits(:,2) == molecule, :);
        
        goodFitFlag(:,molecule) = moleculeFitParam(:,13);
        goodFit = goodFitFlag(:,molecule) > 0;
        goodFitX = goodFit;
        goodFitY = goodFit;
        goodFitZ = goodFit;
        
        % compute deviation with respect to bead location averaged over last
        % numSyncFrames frames of the movie
        devX(:,molecule) = moleculeFitParam(:,21) ...
            - mean(moleculeFitParam(any(bsxfun(@eq,moleculeFitParam(:,1), syncFrames),2),21));
        devY(:,molecule) = moleculeFitParam(:,22) ...
            - mean(moleculeFitParam(any(bsxfun(@eq,moleculeFitParam(:,1), syncFrames),2),22));
        devZ(:,molecule) = moleculeFitParam(:,23) ...
            - mean(moleculeFitParam(any(bsxfun(@eq,moleculeFitParam(:,1), syncFrames),2),23));
        numPhotons(:,molecule) = moleculeFitParam(:,17);
        
        %         devX(:,molecule) = moleculeFitParam(:,21);
        %         devY(:,molecule) = moleculeFitParam(:,22);
        %         devZ(:,molecule) = moleculeFitParam(:,23);
        
        % Filter out activation Frames
        %     goodFitX = abs(devX)<400;
        %     goodFitY = abs(devY)<400;
        %     goodFitZ = abs(devZ)<200;
        
        hDevs = figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
        subplot(2,1,1);
        plot(moleculeFitParam((goodFit & goodFitX),1),devX((goodFit & goodFitX),molecule),'r');
        hold on;
        plot(moleculeFitParam((goodFit & goodFitY),1),devY((goodFit & goodFitY),molecule),'b');
        axis tight;
        legend('x','y');
        title(['Fiduciary ' num2str(molecule)]);
        xlabel('Frame #');
        ylabel('Position (nm)');
        
        subplot(2,1,2);
        plot(moleculeFitParam((goodFit & goodFitZ),1),devZ((goodFit & goodFitZ),molecule),'k');
        axis tight;
        legend('z');
        title(['Fiduciary ' num2str(molecule)]);
        xlabel('Frame #');
        ylabel('Position (nm)');
        print(hDevs,'-djpeg',[outputFilePrefix{stack} 'fiduciary ' num2str(molecule) ' stats.jpg']);
        
        % write fiduciary data to Excel spreadsheet
%         xlswrite([outputFilePrefix{stack} 'fiduciary deviations.xlsx'], ...
%             [textHeader; num2cell([(1:numFrames)' devX(:,molecule) devY(:,molecule) ...
%             devZ(:,molecule) goodFitFlag(:,molecule) numPhotons(:,molecule)])], ...
%             ['fiduciary ' num2str(molecule)]);
        
        % if particle was fit successfully, add its movement to the average
        avgDevX = avgDevX + (goodFit & goodFitX).*devX(:,molecule);
        avgDevY = avgDevY + (goodFit & goodFitY).*devY(:,molecule);
        avgDevZ = avgDevZ + (goodFit & goodFitZ).*devZ(:,molecule);
        numValidFits = numValidFits + (goodFit & goodFitX & goodFitY & goodFitZ);
    end
    
    %% compute average movement of all fiduciaries
    avgDevX = avgDevX./numValidFits;
    avgDevY = avgDevY./numValidFits;
    avgDevZ = avgDevZ./numValidFits;
    
%     figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
%     subplot(2,1,1);
%     plot(avgDevX,'r');
%     hold on;
%     plot(avgDevY,'b');
%     axis tight;
%     legend('x','y');
%     title('Average of all fiduciaries');
%     xlabel('Frame #');
%     ylabel('Position (nm)');
%     
%     subplot(2,1,2);
%     plot(avgDevZ,'k');
%     axis tight;
%     legend('z');
%     title('Average of all fiduciaries');
%     xlabel('Frame #');
%     ylabel('Position (nm)');
%     print('-djpeg',[outputFilePrefix{stack} 'average of all fid stats.jpg']);
    
    % output fiduciary correction data
    save([outputFilePrefix{stack} 'fiduciary corrections.mat'],'avgDevX','avgDevY',...
        'avgDevZ','numValidFits','devX','devY','devZ','goodFitFlag','numPhotons');
    
    % output fiduciary correction data to Excel spreadsheet
%     textHeader = {'frame number' 'deviation in x (nm)' 'deviation in y (nm)' ...
%         'deviation in z (nm)' 'number of valid fiduciaries in average'};
%     xlswrite([outputFilePrefix{stack} 'fiduciary deviations.xlsx'], ...
%         [textHeader; num2cell([(1:numFrames)' avgDevX avgDevY avgDevZ ...
%         numValidFits])], 'average fiduciary deviations');
    
    
    
end


end
