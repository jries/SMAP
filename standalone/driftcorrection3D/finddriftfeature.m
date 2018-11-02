function [drift,driftinfo]=finddriftfeature(pos,par)
global SMAP_stopnow
%pos.xnm .ynm .frame
%frame starts with 1, ascending order, maybe not necessary?
% XXX locData: sort function, general
% parameters and typical values (please adopt)
% par.drift_pixrec=15; %pixelsize of reconstructed images in nm
% par.drift_window=11; %size of region in pixels which gets fittd to determine
% displacement
% par.drift_timepoints=10; %number of time points evaluated 
% par.drift_maxdrift=500; %maximal drift in nm (not crucial, rather choose to
% high
%par.drift_maxpix maximum size of recontsructed image

%other functions needed:
%myhist2
%my2DGaussfit

%copyright: Jonas Ries, EMBL, jonas.ries@embl.de
if isfield(par,'repetitionname')
    rn=par.repetitionname;
else
    rn='';
end
results_ax1=initaxis(par.resultstabgroup,['CC' rn]);
results_ax3=initaxis(par.resultstabgroup,['normalized CC' rn]);

%here: rather from par, from channel range. Make sure it does not get
%displaced. Fill in outside.


lastframe=round(par.framestop);
firstframe=round(par.framestart);
numframes=lastframe-firstframe+1;

%% calculate movie and FFT of movie
pixrec=par.drift_pixrec; %in nm
window=ceil((par.drift_window-1)/2);
timepoints=par.drift_timepoints; %how many timepoints
maxdrift=par.drift_maxdrift; %in nanometers

mx=[min(pos.xnm) max(pos.xnm)]; %ROI which is used for drift correction. 
my=[min(pos.ynm) max(pos.ynm)]; %You can put your own routine here

srec(1)=round((mx(2)-mx(1))/pixrec);
srec(2)=round((my(2)-my(1))/pixrec);

if max(srec)>par.drift_maxpixels %too large for reconstruction, fold back
    pos.xnm=pos.xnm-min(pos.xnm);pos.ynm=pos.ynm-min(pos.xnm);
    maxnm=par.drift_maxpixels*par.drift_pixrec;
    pos.xnm=mod(pos.xnm,maxnm);pos.ynm=mod(pos.ynm,maxnm);
    mx=[min(pos.xnm) max(pos.xnm)]; %ROI which is used for drift correction. 
    my=[min(pos.ynm) max(pos.ynm)]; %You can put your own routine here

    srec(1)=round((mx(2)-mx(1))/pixrec);
    srec(2)=round((my(2)-my(1))/pixrec);
end

% srim= histrender(posr,mx, my, pixrec, pixrec);

nfftexp=2^ceil(log2(max(max(srec),256))); %for fft use power of 2
if nfftexp>2500
    nfftexp=round(max(srec(1:2))/2)*2;
end
noff=nfftexp/2+1; 
disp('make movie')
Fmovier=makemovie;  %calculate fourier transforms of reconstructed images
disp('find displacement')
[ddx, ddy,errx,erry]= finddisplacements2; %determine displacements

ddx=ddx*pixrec; %convert displacements into nm
ddy=ddy*pixrec;

%% bin displacements
[dx,sdxc]=bindisplacementfit(ddx,errx); %determine displacement for each time point
[dy,sdyc]=bindisplacementfit(ddy,erry);

% dx0h=[0; dx];
% dy0h=[0; dy];
s=size(ddx);
% ddxplot=cumsum(diff(ddx));
% ddyplot=cumsum(diff(ddy));
ddxplot=ddx;
ddyplot=ddy;
for kn=1:s(1)
    ddxplot(:,kn)=ddx(:,kn)-ddx(kn,kn)+dx(kn);
    ddyplot(:,kn)=ddy(:,kn)-ddy(kn,kn)+dy(kn);
end




%interpolate displacemnt for all frames
cfit1=(0:length(dx)-1)*binframes+binframes/2+firstframe; %positions of time points
ctrue=(1:par.maxframeall)'; %positions of frames
if length(dx)>9
    [~,sdx,inlier,outlier]=robustMean(ddxplot,2,15);%std for each time point, used for interpolation
    [~,sdy]=robustMean(ddyplot,2,15);
    sdxm=robustMean(sdx)/2;
    sdx(sdx<sdxm)=sdxm;
    sdym=robustMean(sdy)/2;
    sdy(sdy<sdym)=sdym;
    sdxm=robustMean(sdx);
    sdym=robustMean(sdy);
%     indgx=sdx<5*sdxm;
%     indgy=sdy<5*sdym;
end
if length(dx)<=9
    sdx=std(ddxplot,0,2); %std for each time point, used for interpolation
    sdy=std(ddyplot,1,2);
%     indgx=true(size(dx));
%     indgy=true(size(dy));
end
 indgx=true(size(dx));
    indgy=true(size(dy));



% sdx(sdx>5*sdxm)=inf;
% sdy(sdy>5*sdym)=inf;


wx=1./sdx.^2;
wy=1./sdy.^2;


%add to weights
% wx=wx+0.0*mean(wx);
% wy=wy+0.0*mean(wy);

% wx(end)=wx(end)/4;
% wy(end)=wy(end)/4;

% h=cfit1(2)-cfit1(1)
% pset=1/(1+h^3/6)

%give higher weight to first data point:
% wx(1)=wx(1)*2;
% wy(1)=wy(1)*2;
 h=cfit1(2)-cfit1(1);
 if ~isempty(par.smoothpar)
    pset=1/(1+h^3/60*par.smoothpar);
 else
     pset=[];
 end
switch par.smoothmode.Value
    case 1 % smoothing spline
        [dxt,px] = csaps(double(cfit1(indgx)),double(dx(indgx)),double(pset),double(ctrue),wx(indgx)) ;
        [dyt,py] = csaps(double(cfit1(indgy)),double(dy(indgy)),double(pset),double(ctrue),wy(indgy)) ;
    
    case 2
        dxt = interp1(double(cfit1(indgx)),double(dx(indgx)),double(ctrue)) ;
        dyt = interp1(double(cfit1(indgy)),double(dy(indgy)),double(ctrue)) ;
end
framesall=(1:par.maxframeall);%-firstframe+1;
binend=floor(1*binframes/2);
% dxtt=zeros((par.maxframeall),1);dytt=dxtt;
dxtt=dxt;
dxtt(1:firstframe-1+binframes/2)=dxtt(firstframe-1+binframes/2+1);
% dxtt(firstframe:lastframe)=dxt;
dxtt(lastframe+1-binend:end)=dxtt(lastframe+1-binend);

dytt=dyt;
dytt(1:firstframe-1+binframes/2)=dytt(firstframe-1+binframes/2+1);
% dytt(firstframe:lastframe)=dyt;
dytt(lastframe+1-binend:end)=dytt(lastframe+1-binend);


results_ax2=initaxis(par.resultstabgroup,['dxy/frame' rn]);

subplot(1,2,1)
hold off
plot(ddxplot)
hold on
plot(dx,'k','LineWidth',1.5);
plot(sdx,'k:')
sx=(max(dx)-min(dx));
ylim([min(dx)-sx/2 max(dx)+sx/2])
axis tight

subplot(1,2,2)
hold off
plot(ddyplot)
hold on
plot(dy,'k','LineWidth',1.5);
plot(sdy,'k:')
sy=(max(dx)-min(dx));
ylim([min(dy)-sy/2 max(dy)+sy/2])
axis tight

if par.drift_reference
    dxtt=dxtt-dx(end-1);
    dytt=dytt-dy(end-1);
    dx=dx-dx(end-1);
    dy=dy-dy(end-1);
end
    

driftinfo.dx=dx;
driftinfo.dy=dy;
driftinfo.dxplot=ddxplot;
driftinfo.dyplot=ddyplot;
driftinfo.dxt=dxtt;
driftinfo.dyt=dytt;
driftinfo.binframes=cfit1;
%
initaxis(par.resultstabgroup,['dxy/frame final' rn]);

hold off
plot(cfit1,dx,'x',framesall,dxtt,'k')
hold on
plot(cfit1,dy,'o',framesall,dytt,'r')
xlabel('frame')
ylabel('dx, dy (nm)')
drawnow

initaxis(par.resultstabgroup,['dx vs dy' rn]);
hold off
plot(dxtt,dytt,'k')
hold on
plot(dx,dy,'ro')
plot(dx(1),dy(1),'gx')
xlabel('dx')
ylabel('dy')
drawnow
axis equal

drift.x=dxtt;
drift.y=dytt;

% asdafd
% fitposc=adddrift(positions,dxt,dyt); %recalculate positions

function Fmovier=makemovie %calculate fourier transforms of images
%     posr.x=pos.xnm;posr.y=pos.ynm;
    binframes=2*ceil(numframes/timepoints/2+1);
    frameranges=[firstframe:binframes:lastframe lastframe] ;  
    timepoints=length(frameranges)-1;
    Fmovier=zeros(nfftexp,nfftexp,timepoints,'single');
    for k=1:timepoints
        indframe=pos.frame<frameranges(k+1)&pos.frame>=frameranges(k);
        posr.x=pos.xnm(indframe);posr.y=pos.ynm(indframe);
        imager=histrender(posr,mx, my, pixrec, pixrec)';
        Fmovier(:,:,k)=fft2(imager,nfftexp,nfftexp);
        if SMAP_stopnow
            error('execution stopped by user');
        end
%         figure(89)
%         imagesc(imager)
%         waitforbuttonpress
    end
end

function [ddx, ddy,errx,erry]= finddisplacements2 % find displacements
s=size(Fmovier);
dnumframesh =s(3);
ddx=zeros(dnumframesh-1);ddy=zeros(dnumframesh-1);
errx=ddx;
erry=ddy;
% fhold=imagesc(1,'Parent',results_ax1);
timerh=tic;
for k=1:dnumframesh-1
    
    for l=k+1:dnumframesh
        cc=Fmovier(:,:,k).*conj(Fmovier(:,:,l));
        ccf=fftshift(ifft2(cc));
        [mx,my,outim,outimnorm,errx(k,l),erry(k,l)]=findmaximumgauss(real(ccf),window); %maximum by Gaussian fitting
        dxh=mx-noff; dyh=my-noff;
        ddx(k,l)=dxh; ddy(k,l)=dyh;
        ddx(l,k)=-dxh; ddy(l,k)=-dyh;
        errx(l,k)=errx(k,l);erry(l,k)=erry(k,l);
    
    if isfield(par,'showresults') && par.showresults && toc(timerh)>0.5
        timerh=tic;
        fhold=imagesc(outim,'Parent',results_ax1);
        imagesc(outimnorm,'Parent',results_ax3)
        results_ax3.Title.String=num2str(k/dnumframesh+(l-k)/dnumframesh^2);
        results_ax1.Title.String=num2str(k/dnumframesh+(l-k)/dnumframesh^2);
        drawnow
        if SMAP_stopnow
            error('execution stopped by user');
        end
    end
    end
%     disp(k/dnumframesh)
end
end

function [x,y,outim,outimnorm,errx,erry]=findmaximumgauss(img,window)
s=size(img);
win=maxdrift/pixrec; %maxdrift
cent=round(max(1,s(1)/2-win):min(s(1)/2+win,s(1)));
imfm=img(cent,cent);
imfm=filter2(ones(5)/5^2,imfm); %filter a little for better maximum search
[inten,ind]=max(imfm(:)); %determine pixel with maximum intensity to center roi for fitting
[mxh,myh]=ind2sub(size(imfm),ind);
mxh=mxh+cent(1)-1;
myh=myh+cent(1)-1;
%now determine maximum
smallframe=double(img(mxh-window:mxh+window,myh-window:myh+window));
[fitout,outim,outimnorm,ci]=my2Dgaussfit(smallframe,[window+1,window+1,inten,min(smallframe(:)),max(2,3/window),max(2,3/window),0],3);
x=mxh-window+fitout(1)-1;y=myh-window+fitout(2)-1;
dc=ci(:,2)-ci(:,1);
errx=dc(1);erry=dc(2);
end
end

function [dx2,sdx2]=bindisplacementfit(ddx,errx)
% sf=size(ddx);
%idea: we measure displacements between every frames (dxik=xi-xk, xi is 
%displacement for frame i). Use all xi as fit parameters, fit function
%calculates dxik. Robust fit.
% weights=1./(errx+.1).^2;
%startp
difddx=diff(ddx);
fp0=cumsum(median(difddx,2));
% fp0=zeros(sf(1)-1,1);
options=statset('nlinfit');
options=statset(options,'Robust','on');
[fp,r,J,COVB,mse] = nlinfit(ddx(:),ddx(:),@bindispf,fp0,options);
ci = nlparci(fp,r,'covar',COVB);
 dx2=[0; fp];
%  dx2=fp;
sx=ci(:,2)-ci(:,1);
sdx2=[mean(sx);sx];
end

function out=bindisp(fp,ddx)
fph=[0; fp];
ddxf=zeros(length(fph));
for k=1:length(fph)
    ddxf(k,:)=fph(k)-fph; %calculate difference
end


% out=ddxf(:);
out=ddxf;
end
function out=bindispf(fp,ddx)
out=bindisp(fp);
out=out(:);
end
