function plotTimeProfiles4Pi
global pathh plotfig plotfft plotGaussfit

addSMAPpath;
[f,pathh]=uigetfile([pathh '*.dcimg']);
if ~f
    return
end

plotfig=[];plotfft=[];plotGaussfit=[];

il=imageloaderAll([pathh f]);
images=il.getmanyimages([],'mat');
images=images-myquantilefast(images(:),0.05,1000000);
numf=size(images,3);
fg=figure(127);
fg.Position(3)=1700; fg.Position(1)=10;
ax=axes(fg);
ax.Position=[0.05, 0.15,0.9,.8];


h.ax=ax;
h.findmax=uicontrol('Style','checkbox','String','center max','Units','normalized','Position',[0.825,0.05,0.05,0.05],'Value',1,'Parent',fg);
h.slider=uicontrol('Style','slider','Units','normalized','Position',[0.05,0.05,0.5,0.05],'Max',numf,'Value',1,'Min',1,'Parent',fg);
h.framenum=uicontrol('Style','edit','Units','normalized','Position',[0.55,0.05,0.025,0.05],'String','1','Parent',fg);
h.fitmodel=uicontrol('Style','popupmenu','String',{'none','SxSy','Sx'},'Units','normalized','Position',[0.6,0.05,0.05,0.05],'Parent',fg);


h.roisize=uicontrol('Style','edit','String','5','Units','normalized','Position',[0.875,0.05,0.025,0.05],'Parent',fg);


h.plotfft=uicontrol('Style','checkbox','Units','normalized','Position',[0.675,0.05,0.05,0.05],'String','FFT','Value',0,'Parent',fg);
h.fftrate=uicontrol('Style','edit','Units','normalized','Position',[0.7,0.05,0.025,0.05],'String','20','Parent',fg);
h.globalc=uicontrol('Style','checkbox','Units','normalized','Position',[0.725,0.05,0.05,0.05],'String','Glob contr','Parent',fg);
h.trendline=uicontrol('Style','checkbox','Units','normalized','Position',[0.775,0.05,0.05,0.05],'String','trend','Parent',fg);
h.trendlinep=uicontrol('Style','edit','Units','normalized','Position',[0.805,0.05,0.02,0.05],'String','1','Parent',fg);

h.plot=uicontrol('Style','pushbutton','String','plot','Units','normalized','Position',[0.9,0.05,0.05,0.05],'Callback',{@plot_callback,images,h},'Parent',fg);
h.slider.Callback={@slider_callback,images,h};
h.framenum.Callback={@slider_callback,images,h};
% plotfig=figure;
slider_callback(h.slider,0,images,h)
end

function slider_callback(a,b,images,h)
if strcmp(a.Style,'slider')
    fr=ceil(a.Value);
else
    fr=ceil(str2double(a.String));
end
implot=images(:,:,fr);
if h.globalc.Value
    implot(1,1)=max(images(:));
end
imagesc(h.ax,implot)
axis(h.ax,'equal')
h.framenum.String=num2str(fr);
h.slider.Value=fr;
end

function plot_callback(a,b,images,h)
global plotfig  plotfft plotGaussfit
sigma=0.5;
pp=round(h.ax.CurrentPoint(1,1:2));
fr=ceil(h.slider.Value);
if h.findmax.Value
    hs=fspecial('gaussian',5,sigma);
    region=filter2(hs,images(pp(2)-5:pp(2)+5,pp(1)-5:pp(1)+5,fr));
    [~,indm]=max(region(:));
    [x,y]=ind2sub(size(region),indm);
    pp(2)=pp(2)-6+x;
    pp(1)=pp(1)-6+y;
end
if isempty(plotfig) || ~isvalid(plotfig)
    plotfig=figure;
end
f=figure(plotfig);
rois=round(str2double(h.roisize.String));
dr=round(rois/2);
r=-dr:-dr+rois;
intensity=squeeze(sum(sum(images(pp(2)+r,pp(1)+r,:),1),2));
 x=(1:length(intensity))';
plot(x,intensity,'-')
hold on
plot(0,0,'.')
xlabel('frame')
ylabel('intensity')
if h.trendline.Value
    [ptrndl,p]=csaps(x,intensity);
    
    pint=p*str2double(h.trendlinep.String);
    fp=fit(x,intensity,'smoothingspline','SmoothingParam',pint);
    plot(x,fp(x));
    
end


if h.plotfft.Value
    if isempty(plotfft) || ~isvalid(plotfft)
        plotfft=figure; hold on
    end
   
    pwSpec=figure(plotfft);
    Fs = str2double(h.fftrate.String);
    t = 0:1/Fs:(length(intensity)-1)*(1/Fs);

    N = length(intensity);

    xdft = fft(intensity);
    xdft = xdft(1:N/2+1);
    psdx = (1/(Fs*N)) * abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);
    freq = 0:Fs/length(intensity):Fs/2;
    [~,freqb]=bintrace(freq,10);
    [~,psdxb]=bintrace(psdx,10);
    plot(freqb,10*log10(psdxb))
    grid on
    title('Periodogram Using FFT')
    xlabel('Frequency (Hz)')
    ylabel('Power/Frequency (dB/Hz)')
end

%Gaussian fit:
switch h.fitmodel.Value
    case 1 % none
    case 2 %SXSY
            if isempty(plotGaussfit) || ~isvalid(plotGaussfit)
                plotGaussfit=figure;
            end
        roisize2=4; %effecgive 9
        r=-roisize2:roisize2;
        imagesfit=single(images(pp(2)+r,pp(1)+r,:));
        P=mleFit_LM(imagesfit,4,30,1.5);
        sx=P(:,5);
        sy=P(:,6);
        figure(plotGaussfit);
        subplot(3,2,1);
        w=sx.^2-sy.^2;
        plot(x,w)
        xlabel('frame')
        ylabel('sx^2-sy^2');
        hold on
        
        %% FFT sigma
%         pwSpecsigma=figure;
    Fs = str2double(h.fftrate.String);
    t = 0:1/Fs:(length(sx)-1)*(1/Fs);

    N = length(sx);

    xdftw = fft(w);
    
    
      P2w = abs(xdftw/N);
P1w = P2w(1:N/2+1);
P1w(2:end-1) = 2*P1w(2:end-1);
    f = Fs*(0:(N/2))/N;
    subplot(3,2,2)
plot(f(2:end),P1w(2:end)) 
title('Single-Sided Amplitude Spectrum of W(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
    
    
    
%     xdftw = xdftw(1:N/2+1);
%     psdxw = (1/(Fs*N)) * abs(xdftw).^2;
%     psdxw(2:end-1) = 2*psdxw(2:end-1);
%     freq = 0:Fs/length(sx):Fs/2;
%     [~,freqb]=bintrace(freq,10);
%     [~,psdxbw]=bintrace(psdxw,10);
%     plot(freqb,10*log10(psdxbw))
%     grid on
%     title('w Using FFT')
%     xlabel('Frequency (Hz)')
%     ylabel('Power/Frequency (dB/Hz)')
        
%     %% power vx
%     vx = (P(2:end,1)-P(1:end-1,1))*128/0.003;
%     subplot(4,2,3)
%         plot(x(1:end-1),vx)  
%          xlabel('frame')
%         ylabel('vx');
%         hold on
%     
%% plot xy
%     figure;
        subplot(3,2,3)
        plot(x,P(:,1)-mean(P(:,1)),'b') 
        hold on
        plot(x,P(:,2)-mean(P(:,2)),'r')
         xlabel('frame')
        ylabel('x blue, y red');
        hold on
    
%     pwSpecvx=figure;
    Fs = str2double(h.fftrate.String);
    t = 0:1/Fs:(length(sx)-2)*(1/Fs);

    N = length(sx);

 %% plot fft xy
    xdftvx = (fft(P(:,1)));
%     
%      dF = Fs/N;                      % hertz
%    f = -Fs/2:dF:Fs/2-dF;           % hertz
%    figure;
%    plot(f,abs(xdftvx)/N);
%    xlabel('Frequency (in hertz)');
%    ylabel('vx)');
%    title('Magnitude Response');
%     
    P2 = abs(xdftvx/N);
P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(N/2))/N;
    subplot(3,2,4)                              
plot(f(2:end),P1(2:end),'b') 
title('Single-Sided Amplitude Spectrum of XY(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
    hold on

xdftvx = (fft(P(:,2)));
%     
%      dF = Fs/N;                      % hertz
%    f = -Fs/2:dF:Fs/2-dF;           % hertz
%    figure;
%    plot(f,abs(xdftvx)/N);
%    xlabel('Frequency (in hertz)');
%    ylabel('vx)');
%    title('Magnitude Response');
%     
    P2 = abs(xdftvx/N);
P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(N/2))/N;
    subplot(3,2,4)  
    hold on
plot(f(2:end),P1(2:end),'r') 

% %% tilt xy 45degree

%% plot xy tilt 45 degree
%     figure;
        subplot(3,2,5)
        xt = (P(:,1)-P(:,2))*sin(pi/4);
        plot(x,xt-mean(xt),'b') 
        hold on
        yt = (P(:,1)+P(:,2))*sin(pi/4)
        plot(x,yt-mean(yt),'r') 
         xlabel('frame')
        ylabel('tilt x blue, y red');
        hold on
    
%     pwSpecvx=figure;
    Fs = str2double(h.fftrate.String);
    t = 0:1/Fs:(length(sx)-2)*(1/Fs);

    N = length(sx);

xdftvx = (fft(xt));
%     
%      dF = Fs/N;                      % hertz
%    f = -Fs/2:dF:Fs/2-dF;           % hertz
%    figure;
%    plot(f,abs(xdftvx)/N);
%    xlabel('Frequency (in hertz)');
%    ylabel('vx)');
%    title('Magnitude Response');
%     
    P2 = abs(xdftvx/N);
P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(N/2))/N;
    subplot(3,2,6)                              
plot(f(2:end),P1(2:end),'b') 
title('Single-Sided Amplitude Spectrum of XYt(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
    hold on

xdftvx = (fft(yt));
%     
%      dF = Fs/N;                      % hertz
%    f = -Fs/2:dF:Fs/2-dF;           % hertz
%    figure;
%    plot(f,abs(xdftvx)/N);
%    xlabel('Frequency (in hertz)');
%    ylabel('vx)');
%    title('Magnitude Response');
%     
    P2 = abs(xdftvx/N);
P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(N/2))/N;
    subplot(3,2,6)  
    hold on
plot(f(2:end),P1(2:end),'r') 





% title('Single-Sided Amplitude Spectrum of Y(t)')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')

%     xdftvx = xdftvx(1:N/2+1);
%     psdxvx = (1/(Fs*N)) * abs(xdftvx);
%     psdxvx(2:end-1) = 2*psdxvx(2:end-1);
% % psdxvx = xdftvx;
%     freq = 0:Fs/(length(sx)-1):Fs/2;
%     [~,freqb]=bintrace(freq,10);
%     [~,psdxbxv]=bintrace(psdxvx,10);
% %     plot(freqb,10*log10(psdxbxv))
%     plot(freqb,psdxbxv)
%     grid on
%     title('vx Using FFT')
%     xlabel('Frequency (Hz)')
%     ylabel('amplitituted')
        %%
%         figure(plotGaussfit);
%         w=sx;
%         plot(x,w)  
%          xlabel('frame')
%         ylabel('sx');
%         hold on
%         
        
        
        
        
    case 3 %Sx only
        if isempty(plotGaussfit) || ~isvalid(plotGaussfit)
            plotGaussfit=figure;
        end
        roisize2=4; 
        r=-roisize2:roisize2;
        imagesfit=single(images(pp(2)+r,pp(1)+r,:));
        P=mleFit_LM(imagesfit,2,30,1.5);
        sx=P(:,5);
     
        
        
        
        
end
% ax=gca;
% ax.YLim(1)=0;
end