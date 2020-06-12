% [kcMax,A0,kcGM,d0,d] = getDcorr(im,r,Ng,figID)
% ---------------------------------------
%
% Estimate the image cut-off frequency based on decorrelation analysis
%
% Inputs:
%  im        	2D image to be analyzed
%  r           	Fourier space sampling of the analysis (default: r = linspace(0,1,50)
%  Ng			Number of high-pass filtering (default: Ng = 10)
%  figID		If figID > 1, curves will be plotted in figure(figID)
%               if figID == 'fast', enable fast resolution estimate mode
%
% Outputs:
%  kcMax        Estimated cut-off frequency of the image in normalized frequency
%  A0			Amplitude of the local maxima of d0
%  kcGM			Estimated cut-off frequency using Geometric-Mean metric
%  d0 			Decorrelation function before high-pass filtering
%  d			All decorrelation functions
%
% ---------------------------------------
%
% A detailled description of the method can be found in : 
% "Descloux, A., K. S. Grußmayer, and A. Radenovic. "Parameter-free image 
% resolution estimation based on decorrelation analysis."
% Nature methods (2019): 1-7."
% 
%   Copyright © 2018 Adrien Descloux - adrien.descloux@epfl.ch,
%   Ecole Polytechnique Federale de Lausanne, LBEN,
%   BM 5.134, Station 17, 1015 Lausanne, Switzerland.
%
%  	This program is free software: you can redistribute it and/or modify
%  	it under the terms of the GNU General Public License as published by
% 	the Free Software Foundation, either version 3 of the License, or
%  	(at your option) any later version.
%
%  	This program is distributed in the hope that it will be useful,
%  	but WITHOUT ANY WARRANTY; without even the implied warranty of
%  	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  	GNU General Public License for more details.
%
% 	You should have received a copy of the GNU General Public License
%  	along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [kcMax,A0,kcGM,d0,d] = getDcorr(im,r,Ng,figID)

if nargin < 4; figID = 0; end
if ischar(figID)
    figID = 0;
    fastMode = 1;
else
    fastMode = 0;
end
if nargin < 3; Ng = 10;  end
if nargin < 2; r = linspace(0,1,50); end

% input check
if length(r) < 30
    r = linspace(min(r),max(r),30);
end
if Ng < 5
    Ng = 5;
end
%%
im = single(im);
im = im(1:end-not(mod(size(im,1),2)),1:end-not(mod(size(im,2),2))); % odd number of pixels

[X,Y] = meshgrid(linspace(-1,1,size(im,2)),linspace(-1,1,size(im,1)));
R = sqrt(X.^2 + Y.^2);
Nr = length(r);
if isa(im,'gpuArray')
    r = gpuArray(r); 
    R = gpuArray(R);
end
% In : Fourier normalized image
In = fftshift(fftn(fftshift(im))); In = In./abs(In); In(isinf(In)) = 0; In(isnan(In)) = 0;
mask0 = R.^2 < 1^2;
In = mask0.*In; % restric all the analysis to the region r < 1

if figID
    fprintf('Computing dcorr: ');
end
% Ik : Fourier transform of im
Ik = mask0.*fftshift(fftn(fftshift(im)));

c = sqrt(sum(sum(Ik.*conj(Ik))));

r0 = linspace(r(1),r(end),Nr);

for k = length(r0):-1:1
    cc = getCorrcoef(Ik,(R.^2 < r0(k)^2).*In,c);
    if isnan(cc); cc = 0; end
    d0(k) = gather(cc); % gather if input image is gpuArray 
    if fastMode == 1
        [ind0,snr0] = getDcorrLocalMax(d0(k:end));
        ind0 = ind0 + k-1;
        if ind0 > k % found a local maxima, skip the calculation
            break;
        end
    end
end
if fastMode == 0
    ind0 = getDcorrLocalMax(d0(k:end));
    snr0 = d0(ind0);
end
res0 = gather(r(ind0));

gMax = 2/r0(ind0);
if isinf(gMax); gMax = max(size(im,1),size(im,2))/2;end

% search of highest frequency peak
g = [size(im,1)/4, exp(linspace(log(gMax),log(0.15),Ng))];
d = zeros(Nr,2*Ng); kc = res0; SNR = snr0; gm = res0;
if fastMode == 0
    ind0 = 1;
else
    if ind0 > 1
        ind0 = ind0 -1;
    end
end
for refin = 1:2 % two step refinement

    for h = 1:length(g)
        Ir = Ik.*(1 - exp(-2*g(h)*g(h)*R.^2)); % Fourier Gaussian filtering
        c = sqrt(sum(sum(abs(Ir).^2)));
        
        for k = length(r):-1:ind0

            if isa(im,'gpuArray')
                cc = getCorrcoef(Ir,In.*(R.^2 < r(k)^2),c);
                if isnan(cc); cc = 0; end
                d(k,h + Ng*(refin-1)) = gather(cc);
            else
                mask = (R.^2 < r(k)^2);
                cc = getCorrcoef(Ir(mask),In(mask),c);
                if isnan(cc); cc = 0; end
                d(k,h + Ng*(refin-1)) = cc;
            end
            if fastMode
                [ind,snr] = getDcorrLocalMax(d(k:end,h + Ng*(refin-1)));
                ind = ind +k-1;
                if ind > k % found a local maxima, skip the calculation
                    break;
                end
            end

        end
        if fastMode == 0
            ind = getDcorrLocalMax(d(k:end,h + Ng*(refin-1)));
            snr = d(ind,h + Ng*(refin-1));
            ind = ind +k-1;
        end
        kc(h + Ng*(refin-1)+1) = gather(r(ind));
        SNR(h + Ng*(refin-1)+1) = snr;
        gm(h + Ng*(refin-1)+1) = sqrt(snr*kc(end));
        if figID
        	fprintf('-');
        end
    end

% refining the high-pass threshold and the radius sampling
    if refin == 1

    % high-pass filtering refinement
        indmax = find(kc == max(kc));
        ind1 = indmax(end);
        if ind1 == 1 % peak only without highpass 
            ind1 = 1;
            ind2 = 2;
            g1 = size(im,1);
            g2 = g(1);
        elseif ind1 >= numel(g)
            ind2 = ind1-1;
            ind1 = ind1-2;
            g1 = g(ind1); g2 = g(ind2);
        else
            ind2 = ind1;
            ind1 = ind1-1;
            g1 = g(ind1); g2 = g(ind2);
        end
        g = exp(linspace(log(g1),log(g2),Ng));
        
        % radius sampling refinement

        r1 = kc(indmax(end))-(r(2)-r(1)); r2 = kc(indmax(end))+0.3;
        if r1 < 0 ; r1 = 0; end
        if r2 > 1; r2 = 1; end
        r = linspace(r1,min(r2,r(end)),Nr);
        ind0 = 1;
        r2 = r;
    end
end
if figID
    fprintf(' -- Computation done -- \n');
end

% add d0 results to the analysis (usefull for high noise images)
kc(end+1) = gather(res0);
SNR(end+1) = snr0;

% % need at least 0.05 of SNR to be even considered
kc(SNR < 0.05) = 0;
SNR(SNR < 0.05) = 0;

snr = SNR;

% output results computation
if ~isempty(kc)
    % highest resolution found 
    [kcMax,ind] = max(kc);
    AMax = SNR(ind);

    % compute the geometric mean to determine the best res/SNR curve
    gm = sqrt(kc.*SNR);
    [~,ind] = max(gm);

    kcGM = kc(ind);
    A0 = snr0; % average image contrast has to be estimated from original image
else
    kcMax = r(2);
    Amax = 0;
    res = r(2);
    A0 = 0;
end

% results display if figID specified
if figID 
    radAv = getRadAvg(log(abs(gather(Ik))+1));
    lnwd = 1.5;
    figure(figID);
    plot(r0,d(:,1:Ng),'color',[0.2 0.2 0.2 0.5]);
    hold on
    radAv(1) = radAv(2); %for plot 
    radAv(end) = radAv(end-1);
    plot(linspace(0,1,length(radAv)),linmap(radAv,0,1),'linewidth',lnwd,'color',[1 0 1])
    for n = 1:Ng
        plot(r2,d(:,n+Ng),'color',[0 0 (n-1)/Ng])
    end
    plot(r0,d0,'linewidth',lnwd,'color','g')
    plot([kcMax kcMax],[0 1],'k')
    for k = 1:length(kc)
        plot(kc(k),snr(k),'bx','linewidth',1)
    end
    hold off
    title(['Dcor analysis : res ~ ',num2str(kcMax,4),', SNR ~ ',num2str(A0,4)])
    xlim([0 1]); ylim([0 1])
    xlabel('Normalized spatial frequencies')
    ylabel('C.c. coefficients')
end