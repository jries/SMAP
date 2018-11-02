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

function f_debugMoleculeFits(totalPSFfits)
        [debugFile, debugPath] = uiputfile({'*.csv';'*.txt'},'Specify a path and name for debug histogram and full data output');
        if isequal(debugFile,0)
            return;
        end
        %makes sure output has an extension
        if debugFile(end-3) ~= '.'
            debugFile = [debugFile '.csv'];
        end
        edges = [-1007:-1000,-3,1,inf];
        vector = histc(totalPSFfits(:,17),edges);
        scrsz = get(0,'ScreenSize');
        hErrors=figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');

        bar(vector(1:end-1))
        set(gca,'XTickLabel',{'LS error', 'amp ratio', 'lobe dist',...
                              'sig ratio', 'sig size', 'pks out',...
                              'amp< 0', 'guess out','fit err',...
                              'good fit'})
        set(gca,'FontSize',8);
                          
        print(hErrors,'-dpng',[debugPath debugFile(1:end-4) ' outcomes 1.png']);
        saveas(hErrors,[debugPath debugFile(1:end-4) ' outcomes.png']);
        % open a file for writing
        [fid,message] = fopen([debugPath debugFile], 'w');
        if ~isempty(message)
            error([debugPath debugFile ': ' message]);
            %return;
        end
        % print a title, followed by a blank line
        fprintf(fid, ['frame num,fit flag,SM idx in frame,template x (pix),template y (pix), template idx,' ...
            'match strength,amp1,amp2,peak1 x (pix),peak1 y (pix),' ...
            'peak2 x (pix),peak2 y (pix), sigma1 (pix),sigma2 (pix),mean bkgnd photons,'...
            'fit error,molecule x (nm),molecule y (nm),DHPSF angle,' ...
            'num photons,interlobe distance,amplitude ratio,sigma ratio,x fid-corrected (nm),y fid-corrected (nm), z fid-corrected (nm),'...
            'photons detected,mean background photons\n']);
        fclose(fid);

        dlmwrite([debugPath debugFile],totalPSFfits(:,[1 17 2:16 18:end]),'-append');
        disp('Debug files written successfully.');
end