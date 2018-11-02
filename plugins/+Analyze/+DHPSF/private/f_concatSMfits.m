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

function [totalPSFfits, numFramesInFiles] = ...
    f_concatSMfits(fitFilePrefix,useFidCorrections,fidFilePrefix)
%clear all;
% close all;
numSyncFrames = 25;
useDenoising = 1;

%% load raw fiduciary data
% [fidFile fidPath] = uigetfile({'*.mat';'*.*'},'Open data file #1 with raw fiduciary fits');

if useFidCorrections
    
%     fidFiles = {};
    
    for fileNum=1:length(fidFilePrefix)
%         fidFiles = [fidFiles; {[fidPath fidFile]}];
        % load data
        load([fidFilePrefix{fileNum} 'raw fits.mat'],'PSFfits','numFrames','numMoles');
        
        if fileNum == 1
            tempPSFfits = PSFfits(:,1:23);
            numFramesInFiles = numFrames;
        else
            numFramesInFiles = [numFramesInFiles numFrames];
            PSFfits(:,1) = PSFfits(:,1) + sum(numFramesInFiles(1:fileNum-1));
            tempPSFfits = [tempPSFfits; PSFfits(:,1:23)];
        end
        
%         fileNum = fileNum+1;
%         [fidFile fidPath] = uigetfile({'*.mat';'*.*'},...
%             ['Open data file #' num2str(fileNum) ' with raw fiduciary fits']);
    end
    PSFfits = tempPSFfits;
    numFrames = sum(numFramesInFiles);
    clear tempPSFfits;
    
    
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
        end
        syncFrames(a)=lastGoodFrame;
        lastGoodFrame = lastGoodFrame - 1;
    end
    
    for molecule = 1:numMoles
        %% extract fitting parameters for this molecule
        moleculeFitParam = PSFfits(PSFfits(:,2) == molecule, :);
        
        goodFitFlag(:,molecule) = moleculeFitParam(:,13);
        goodFit = goodFitFlag(:,molecule) > 0;
        
        % compute deviation with respect to bead location averaged over last
        % numSyncFrames frames of the movie
        devX(:,molecule) = moleculeFitParam(:,21) ...
            - mean(moleculeFitParam(any(bsxfun(@eq,moleculeFitParam(:,1), syncFrames),2),21));
        devY(:,molecule) = moleculeFitParam(:,22) ...
            - mean(moleculeFitParam(any(bsxfun(@eq,moleculeFitParam(:,1), syncFrames),2),22));
        devZ(:,molecule) = moleculeFitParam(:,23) ...
            - mean(moleculeFitParam(any(bsxfun(@eq,moleculeFitParam(:,1), syncFrames),2),23));
        numPhotons(:,molecule) = moleculeFitParam(:,17);
        
        % write fiduciary data to Excel spreadsheet
%         xlswrite([saveFilePrefix 'fiduciary deviations.xlsx'], ...
%             [textHeader; num2cell([(1:numFrames)' devX(:,molecule) devY(:,molecule) ...
%             devZ(:,molecule) goodFitFlag(:,molecule) numPhotons(:,molecule)])], ...
%             ['fiduciary ' num2str(molecule)]);
        
        % if particle was fit successfully, add its movement to the average
        avgDevX = avgDevX + goodFit.*devX(:,molecule);
        avgDevY = avgDevY + goodFit.*devY(:,molecule);
        avgDevZ = avgDevZ + goodFit.*devZ(:,molecule);
        numValidFits = numValidFits + goodFit;
    end
    tempAvgDevX = avgDevX./numValidFits;
    tempAvgDevY = avgDevY./numValidFits;
    tempAvgDevZ = avgDevZ./numValidFits;
    
end

if exist('tempAvgDevX','var')
%     textHeader = {'frame number' 'deviation in x (nm)' 'deviation in y (nm)' ...
%         'deviation in z (nm)' 'number of valid fiduciaries'};
%     
%     xlswrite([saveFilePrefix 'fiduciary deviations.xlsx'], ...
%         [textHeader; num2cell([(1:numFrames)' tempAvgDevX tempAvgDevY  ...
%         tempAvgDevZ numValidFits])], 'average of all fiduciaries');
%     
%     xlswrite([saveFilePrefix 'fiduciary deviations.xlsx'], ...
%         [{'fiduciary files' 'Number of frames'}; ...
%         fidFiles num2cell(numFramesInFiles')], ...
%         'concatenation info');

    % load localization data
%     moleFiles = {};
    for c = 1:fileNum % -1
        % load data
        load([fitFilePrefix{c} 'molecule fits.mat'],'totalPSFfits','numFrames');
        
        if c == 1
            tempPSFfits = totalPSFfits(:,1:27);
        else
            totalPSFfits(:,1) = totalPSFfits(:,1) + sum(numFramesInFiles(1:c-1));
            tempPSFfits = [tempPSFfits; totalPSFfits(:,1:27)];
        end
    end
    totalPSFfits = tempPSFfits;
%     numFrames = sum(numFramesInFiles);
    clear tempPSFfits;
    
    % De-noise the fiduciary tracks
    avgDevX = tempAvgDevX;
    avgDevY = tempAvgDevY;
    avgDevZ = tempAvgDevZ;
    
    [avgDevX_denoised,avgDevY_denoised,avgDevZ_denoised] = f_waveletFidTracks(avgDevX,avgDevY,avgDevZ,1);
    
    % apply fiduciary corrections
    if useDenoising == 1
        xFidCorrected = totalPSFfits(:,25) - avgDevX_denoised(totalPSFfits(:,1));
        yFidCorrected = totalPSFfits(:,26) - avgDevY_denoised(totalPSFfits(:,1));
        zFidCorrected = totalPSFfits(:,27) - avgDevZ_denoised(totalPSFfits(:,1));
    else
        xFidCorrected = totalPSFfits(:,25) - avgDevX(totalPSFfits(:,1));
        yFidCorrected = totalPSFfits(:,26) - avgDevY(totalPSFfits(:,1));
        zFidCorrected = totalPSFfits(:,27) - avgDevZ(totalPSFfits(:,1));
    end
    totalPSFfits = [totalPSFfits xFidCorrected yFidCorrected zFidCorrected];
    clear tempAvgDevX tempAvgDevY tempAvgDevZ;
    
    % output corrected data
%     save([saveFilePrefix 'molecule fits.mat']);
    
    % output excel spreadsheet
%     textHeader = [textHeader {'fiduciary corrected x location (nm)' ...
%         'fiduciary corrected y location (nm)' ...
%         'fiduciary corrected z location (nm)'}];
    
%     xlswrite([saveFilePrefix 'molecule fits.xlsx'], [textHeader; ...
%         num2cell(totalPSFfits)], ...
%         'PSF fits');
%     xlswrite([saveFilePrefix 'molecule fits.xlsx'], ...
%         [{'fiduciary files' 'localization files' 'Number of frames'}; ...
%         fidFiles moleFiles num2cell(numFramesInFiles')], ...
%         'concatenation info');
    
else
    
    % load localization data
%     [moleFile molePath] = uigetfile({'*.mat';'*.*'},'Open data file #1 with PSF localizations');
        
%     moleFiles = {};
    
    for fileNum=1:length(fitFilePrefix)
%         moleFiles = [moleFiles; {[molePath moleFile]}];
        clear numFrames;
        load([fitFilePrefix{fileNum} 'molecule fits.mat'],'totalPSFfits','numFrames');
        
        if fileNum == 1
            tempPSFfits = totalPSFfits(:,1:27);
            if ~exist('numFrames','var')
                numFrames = max(totalPSFfits(:,1));
            end
            numFramesInFiles = numFrames;
        else
            if ~exist('numFrames','var')
                numFrames = max(totalPSFfits(:,1));
            end
            numFramesInFiles = [numFramesInFiles numFrames];
            totalPSFfits(:,1) = totalPSFfits(:,1) + sum(numFramesInFiles(1:fileNum-1));
            tempPSFfits = [tempPSFfits; totalPSFfits(:,1:27)];
        end
        
%         fileNum = fileNum+1;
%         [moleFile molePath] = uigetfile({'*.mat';'*.*'},...
%             ['Open data file #' num2str(fileNum) ' with PSF localizations']);
    end
    totalPSFfits = [tempPSFfits nan(size(tempPSFfits,1),3)];
%     numFrames = sum(numFramesInFiles);
    clear tempPSFfits;
    
    
    
    % output concatenated data
%     save([saveFilePrefix 'molecule fits.mat']);
    
    % output excel spreadsheet
%     xlswrite([saveFilePrefix 'molecule fits.xlsx'], [textHeader; ...
%         num2cell(totalPSFfits)], ...
%         'PSF fits');
%     xlswrite([saveFilePrefix 'molecule fits.xlsx'], ...
%         [{'localization files' 'Number of frames'}; ...
%         moleFiles num2cell(numFramesInFiles')], ...
%         'concatenation info');
    
end
end
