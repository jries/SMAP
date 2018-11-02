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

function [xfid_filt,yfid_filt,zfid_filt] = f_waveletFidTracks(xfid,yfid,zfid,FastVersion)
%denoises raw fiducial tracks by using a combined low-pass/compressing
%filter based on wavelet methods
%
%[xfid_filt,yfid_filt,zfid_filt] =
%f_waveletFidTracks(xfid,yfid,zfid,FastVersion) takes input raw tracks xfid,
%yfid, zfid and outputs filtered tracks xfid_filt, yfid_filt, zfid_filt.
%FastVersion must be either 0 or 1. If FastVersion = 0, the program will
%filter xfid, yfid, zfid one at a time, with user inputs for each set of
%wavelet paramaters. If FastVersion = 1, the program will apply the same
%preset wavelet parameters to all tracks.
%
%Version 1              2012.05.04       Mikael Backlund



    %remove NaN entries before decomposition
    xfid_nonan = xfid(~isnan(xfid));
    yfid_nonan = yfid(~isnan(yfid));
    zfid_nonan = zfid(~isnan(zfid));
    greatjobx = 0;
    greatjoby = 0;
    greatjobz = 0;

    if FastVersion == 0;
        %%
        % x trackblock
        while greatjobx == 0;
            figure, plot(xfid_nonan,'r')
            set(gcf, 'Position', get(0,'Screensize')); % Maximize figure
            title('x raw track')

            prompt = {'level of x decomposition',...
                'x wavelet basis','x coefficient threshold'};
            name = 'Input for x track wavelet desomposition';
            defaultans = {'3','db3','50'};

            input = inputdlg(prompt,name,1,defaultans);
            close

            decx = str2double(input(1)); %decomposition level
            basisx = cell2mat(input(2)); %wavelet basis
            thrx = str2double(input(3)); %threshold for coeffiecient absolute value

            %wavelet decomposition
            [Cx,Lx] = wavedec(xfid_nonan,decx,basisx);
            %compress wavelets with coefficients below threshold. Also preserve low
            %frequency components
            Cx_filt = Cx.*(abs(Cx)>thrx);
            Cx_filt(1:Lx(1)+Lx(2)) = Cx(1:Lx(1)+Lx(2));    %reconstruct signal from filtered wavelets
            xfid_filt_nonan = waverec(Cx_filt,Lx,basisx);
            %put NaN entries back where they were
            xfid_filt = xfid;
            xfid_filt(~isnan(xfid)) = xfid_filt_nonan;

            figure, plot(xfid,'r')
            hold on
            plot(xfid_filt,'k')
            set(gcf, 'Position', get(0,'Screensize')); % Maximize figure
            legend('x raw track','x filtered track')
            title('Inspect overlay. Press any key to continue')
            pause

            oopsx = 1;
            while oopsx == 1;
                
                button = questdlg('Do you like what you see?',...
                    'Make an assessment','yes','no','yes');
                close

                if strcmp(button,'yes') == 1;
                    greatjobx = 1;
                    oopsx = 0;
                elseif strcmp(button,'no') == 1;
                    greatjobx = 0 ;
                    oopsx = 0;
                else
                    error('User canceled program');
                end
            end


        end
        %%
        %y track block
       while greatjoby == 0;
            figure, plot(yfid_nonan,'r')
            set(gcf, 'Position', get(0,'Screensize')); % Maximize figure
            title('y raw track')

            prompt = {'level of y decomposition',...
                'y wavelet basis','y coefficient threshold'};
            name = 'Input for y track wavelet desomposition';
            defaultans = {'3','db3','50'};

            input = inputdlg(prompt,name,1,defaultans);
            close

            decy = str2double(input(1)); %decomposition level
            basisy = cell2mat(input(2)); %wavelet basis
            thry = str2double(input(3)); %threshold for coeffiecient absolute value

            %wavelet decomposition
            [Cy,Ly] = wavedec(yfid_nonan,decy,basisy);
            %compress wavelets with coefficients below threshold. Also preserve low
            %frequency components
            Cy_filt = Cy.*(abs(Cy)>thry);
            Cy_filt(1:Ly(1)+Ly(2)) = Cy(1:Ly(1)+Ly(2));    %reconstruct signal from filtered wavelets
            yfid_filt_nonan = waverec(Cy_filt,Ly,basisy);
            %put NaN entries back where they were
            yfid_filt = yfid;
            yfid_filt(~isnan(yfid)) = yfid_filt_nonan;

            figure, plot(yfid,'r')
            hold on
            plot(yfid_filt,'k')
            set(gcf, 'Position', get(0,'Screensize')); % Maximize figure
            legend('y raw track','y filtered track')
            title('Inspect overlay. Press any key to continue')
            pause

            oopsy = 1;
            while oopsy == 1;
                
                button = questdlg('Do you like what you see?',...
                    'Make an assessment','yes','no','yes');
                close

                if strcmp(button,'yes') == 1;
                    greatjoby = 1;
                    oopsy = 0;
                elseif strcmp(button,'no') == 1;
                    greatjoby = 0 ;
                    oopsy = 0;
                else
                    error('User canceled program');
                end
            end
       end

        %%
            %z track block
       while greatjobz == 0;
            figure, plot(zfid_nonan,'r')
            set(gcf, 'Position', get(0,'Screensize')); % Maximize figure
            title('z raw track')

            prompt = {'level of z decomposition',...
                'z wavelet basis','z coefficient threshold'};
            name = 'Input for z track wavelet desomposition';
            defaultans = {'3','db3','50'};

            input = inputdlg(prompt,name,1,defaultans);
            close

            decz = str2double(input(1)); %decomposition level
            basisz = cell2mat(input(2)); %wavelet basis
            thrz = str2double(input(3)); %threshold for coeffiecient absolute value

            %wavelet decomposition
            [Cz,Lz] = wavedec(zfid_nonan,decz,basisz);
            %compress wavelets with coefficients below threshold. Also preserve low
            %frequencz components
            Cz_filt = Cz.*(abs(Cz)>thrz);
            Cz_filt(1:Lz(1)+Lz(2)) = Cz(1:Lz(1)+Lz(2));    %reconstruct signal from filtered wavelets
            zfid_filt_nonan = waverec(Cz_filt,Lz,basisz);
            %put NaN entries back where thez were
            zfid_filt = zfid;
            zfid_filt(~isnan(zfid)) = zfid_filt_nonan;

            figure, plot(zfid,'r')
            hold on
            plot(zfid_filt,'k')
            set(gcf, 'Position', get(0,'Screensize')); % Maximize figure
            legend('z raw track','z filtered track')
            title('Inspect overlay. Press any key to continue')
            pause

            oopsz = 1;
            while oopsz == 1;
                
                button = questdlg('Do you like what you see?',...
                    'Make an assessment','yes','no','yes');
                close

                if strcmp(button,'yes') == 1;
                    greatjobz = 1;
                    oopsz = 0;
                elseif strcmp(button,'no') == 1;
                    greatjobz = 0 ;
                    oopsz = 0;
                else
                    error('User canceled program');
                end
            end
       end
    elseif FastVersion == 1;
        %%
        %fast version: enter one set of paramaters for all tracks

        dec = 6; %decomposition level
        basis = 'db3'; %wavelet basis
        thr = 50; %threshold for coefficient absolute value


        %wavelet decomposition
        [Cx,Lx] = wavedec(xfid_nonan,dec,basis);
        [Cy,Ly] = wavedec(yfid_nonan,dec,basis);
        [Cz,Lz] = wavedec(zfid_nonan,dec,basis);

        %compress wavelets with coefficients below threshold. Also preserve low
        %frequency components
        Cx_filt = Cx.*(abs(Cx)>thr);
        Cx_filt(1:Lx(1)+Lx(2)) = Cx(1:Lx(1)+Lx(2));
        Cy_filt = Cy.*(abs(Cy)>thr);
        Cy_filt(1:Ly(1)+Ly(2)) = Cy(1:Ly(1)+Ly(2));
        Cz_filt = Cz.*(abs(Cz)>thr);
        Cz_filt(1:Lz(1)+Lz(2)) = Cz(1:Lz(1)+Lz(2));

        %reconstruct signal from filtered wavelets
        xfid_filt_nonan = waverec(Cx_filt,Lx,basis);
        yfid_filt_nonan = waverec(Cy_filt,Ly,basis);
        zfid_filt_nonan = waverec(Cz_filt,Lz,basis);

        %put NaN entries back where they were
        xfid_filt = xfid;
        xfid_filt(~isnan(xfid)) = xfid_filt_nonan;
        yfid_filt = yfid;
        yfid_filt(~isnan(yfid)) = yfid_filt_nonan;
        zfid_filt = zfid;
        zfid_filt(~isnan(zfid)) = zfid_filt_nonan;
        
        figure, plot(xfid,'r')
        hold on
        plot(xfid_filt,'k')
        legend('x raw track','x filtered track')
        title('x overlay')
        
        figure, plot(yfid,'r')
        hold on
        plot(yfid_filt,'k')
        legend('y raw track','y filtered track')
        title('y overlay')
        
        figure, plot(zfid,'r')
        hold on
        plot(zfid_filt,'k')
        legend('z raw track','z filtered track')
        title('z overlay')
        
    else
        error('FastVersion must be 1 or 0');
    end
        
    
end
    
    
    
    