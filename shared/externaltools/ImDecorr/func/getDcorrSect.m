% [kcMax,A0,kcGM,d0,d] = getDcorrSect(im,r,Ng,Na,figID)
% ---------------------------------------
%
% Estimate the image sectorial cut-off frequency based on decorrelation analysis
%
% Inputs:
%  im        	2D image to be analyzed
%  r           	Fourier space sampling of the analysis (default: r = linspace(0,1,50)
%  Ng			Number of high-pass filtering (default: Ng = 10)
%  Na 			Number of sectors
%  figID		If figID > 1, curves will be plotted in figure(figID)
%
% Outputs:
%  kcMax        Estimated cut-off frequency of the image in normalized frequency for each sectors
%  A0			Amplitude of the local maxima of d0 for each sectors
%  kcGM			Estimated cut-off frequency using Geometric-Mean metric for each sectors
%  d0 			Decorrelation function before high-pass filtering for each sectors
%  d			All decorrelation functions for each sectors
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
          
function [kcMax,A0,kcGM,d0,d] = getDcorrSect(im,r,Ng,Na,figID)
            
if nargin < 5; figID = 0; end
if nargin < 4; Na = 2; end % default is 2 angles [0,90]
if nargin < 3; Ng = 10;  end
if nargin < 2; r = linspace(0,1,50); end


a = linspace(-pi,pi,2*Na+1);

im = single(im);

if mod(size(im,2),2) == 0
    im = im(1:end-1,1:end-1);
end

if figID 
    hwait = waitbar(0,'Computing 2D dcorr');
end

Nr = length(r);

[X,Y] = meshgrid(linspace(-1,1,size(im,2)),linspace(-1,1,size(im,1)));
R = sqrt(X.^2 + Y.^2);
th = imrotate(atan2(Y,X),rad2deg(pi/(2*Na)),'crop');
mask0 = R.^2 < 1^2;

% clean fourier spare of input image
Ik = mask0.*fftshift(fftn(fftshift(im)));
imr = real(ifftshift(ifftn(ifftshift(Ik))));
mt = zeros(size(mask0));

% compute dcorr 0 and find its maxima for all angles
for na = 1:Na
    
    % define Fourier mask
    maskA = (th >= a(na) & th <= a(na+1)) | (th >= a(Na+na) & th <= a(Na+na+1)) | R < 0.01;
    mt = mt + na.*maskA;
    Ir = mask0.*maskA.*fftshift(fftn(fftshift(imr)));
    radAv(na,:) = getRadAvg(gather(log(abs(Ir)+1)));

    % Fourier space normalization
    I = Ir./abs(Ir);
    I(isinf(I)) = 0; I(isnan(I)) = 0;
    
    Ir = Ir(1:(end-1)/2);
    c = sqrt(sum(sum(abs(Ir).^2)));

    count = 0;
    r0 = linspace(0,1,Nr); 

    for k = length(r0):-1:1
        rt = r0(k);
        mask  = R.^2 < rt.^2;
        
        temp = mask.*I;
        temp = temp(1:(end-1)/2); % remove the mean & optimize speed
    
        cc = gather(getCorrcoef(Ir,temp,c));
    
        nanmap = isnan(cc);
        cc(nanmap) = 0; % remove nan
    
        d0(k,na) = cc; % gather if input image is gpuArray 
        count = count +1;
        if figID 
            waitbar(0.1*count/(Na*Nr),hwait);
        end
    end
end
for na = 1:Na
    [ind(na),snr0(na)] = getDcorrMax(d0(:,na));
end

% Display the sectorial mask
% figure;imagesc(mask0.*mt)

gMax = 1./r0(min(ind));
mapinf = isinf(gMax);
gMax(mapinf) = 2/r0(2);

% automatic search of best geometric mean
g = exp(linspace(log(gMax),log(0.14),Ng));

d = []; kc = []; SNR = []; gm = []; ind =[];
r2 = repmat(r0,[Na 1]);
for refin = 1:2 % two step refinement
    for na = 1:Na
        
        % define proper Fourier mask
        maskA = (th >= a(na) & th <= a(na+1)) | (th >= a(Na+na) & th <= a(Na+na+1));
        
        for h = 1:length(g)
            imt = imr - imgaussfilt(imr,g(h));
            
            Ir = mask0.*maskA.*fftshift(fftn(fftshift(imt)));
            I = Ir./abs(Ir);
            I(isinf(I)) = 0; I(isnan(I)) = 0;
            
            Ir = Ir(1:(end-1)/2); % remove the mean
            c = sqrt(sum(sum(abs(Ir).^2)));
            
            for k = size(r2,2):-1:1
                rt = r2(na,k);
                mask  = R.^2 < rt.^2;
                temp = mask.*I;
                temp = temp(1:(end-1)/2); % remove the mean
                cc = gather(getCorrcoef(Ir,temp,c));
                
                map = isnan(cc); cc(map) = 0;
                d(k,na,h) = cc; % gather if input image is gpuArray
                count = count+1;
                if figID
                    waitbar(0.1 + 0.9*count/(Nr*Na*Ng*2),hwait);
                end
            end
            
        end
        
        for h =1:length(g)
        	[ind(na,h),snr(na,h)] = getDcorrMax(d(:,na,h));
        end
        
        kc(na,:) = r2(na,ind(na,:));
        SNR(na,:) = snr(na,:);
        gm(na,:) = sqrt(snr(na,:).*r2(na,ind(na,:)));
        
    end
            % refining the high-pass threshold and the radius sampling
        if refin == 1
            [~,ind] = max(gm,[],2);
            map1 = ind == 1;
            mapend = ind == length(g);
            ind(map1) = 2;
            ind(mapend) = size(g,1) - 1;
            
            g1 = g(min(ind)-1); g2 = g(max(ind)+1);
            g = exp(linspace(log(g1),log(g2),Ng));
            
            r2 = [];
            for na = 1:Na
                rmin = kc(na,ind(na)); rmax = kc(na,ind(na))+0.2;
                map0 = rmin < 0; rmin(map0) = 0;
                map1 = rmax > 1; rmax(map1) = 1;
                r2(na,:) = linspace(rmin,rmax,50);
            end
            
            dtemp = d;
        end
end

% keep all data
d = cat(3,dtemp,d);

% need at least 0.05 of SNR to be even considered
kc(SNR < 0.05) = 0;
SNR(SNR < 0.05) = 0;
snr = SNR;

kcMax = zeros(Na,1); kcGM = kcMax; A0 = kcMax;
for na = 1:Na
	if sum(kc(na,:)) > 0
        % highest resolution found 
        [kcMax(na),ind] = max(kc(na,:));
        SNRMax(na) = SNR(na,ind);

        % compute the geometric mean to determine the best res/SNR curve
        gm = sqrt(kc(na,:).*SNR(na,:));
        [~,ind] = max(gm);

        kcGM(na) = kc(na,ind);
        A0(na) = snr0(na); % average image contrast has to be estimated from original image
    else
        kcMax(na) = r(2);
        SNRMax(na) = 0;
        kcGM(na) = r(2);
        SNRout(na) = 0;
    end
end

if figID
    waitbar(1,hwait);
    delete(hwait)
end
if figID
    x = linspace(-1,1,size(im,2)); y = linspace(-1,1,size(im,1));
    figure(figID)
    imagesc(x,y,log(abs(fftshift(fftn(fftshift((im)))))+1))
    hold on
    th = linspace(0,2*pi,2*Na+1)';
    
    plot(cos(th).*[kcMax; kcMax; kcMax(1)],sin(th).*[kcMax; kcMax; kcMax(1)],'w--','linewidth',1.5)
    hold off
    
    figure(figID+1)
    nsub = ceil(sqrt(Na));
    for n =1:Na
        subplot(nsub,nsub,n)
    	lnwd = 1.5;
        r0 = linspace(0,1,size(d,1));
        plot(r0,squeeze(d(:,n,1:Ng)),'color',[0.2 0.2 0.2 0.5]);
    	hold on
      	radAv(n,1) = radAv(n,2); %for plot 
        radAv(n,end) = radAv(n,end-1);
      	plot(linspace(0,1,size(radAv,2)),linmap(radAv(n,:),0,1),'linewidth',lnwd,'color',[1 0 1])
      	for k = 1:Ng
         	plot(r2(n,:),squeeze(d(:,n,k+Ng:end)),'color',[0 0 (n-1)/Ng])
        end
       	plot(r0,d0(:,n),'linewidth',lnwd,'color','g')
      	plot([kcMax(n) kcMax(n)],[0 1],'k')
      	for k = 1:size(kc,2)
            plot(kc(n,k),snr(n,k),'bx','linewidth',1)
        end
        hold off
        title(['R: ',num2str(kcMax(n),2),', A: ',num2str(A0(n),2)])
        xlim([0 1]); ylim([0 1])

    end
    
end
