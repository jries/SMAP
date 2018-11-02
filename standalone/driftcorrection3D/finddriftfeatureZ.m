function [drift,driftinfo]=finddriftfeatureZ(pos,par)
global SMAP_stopnow
%pos.xnm .ynm  .znm .frame
%frame starts with 1, ascending order, maybe not necessary?
% XXX locData: sort function, general
% parameters and typical values (please adopt)
% par.drift_pixrec=15; %pixelsize of reconstructed images in nm
% par.drift_window=11; %size of region in pixels which gets fittd to determine
% displacement
% par.drift_timepoints=10; %number of time points evaluated 
% par.drift_maxdrift=500; %maximal drift in nm (not crucial, rather choose to
% high

%other functions needed:
%myhist2
%my2DGaussfit

%copyright: Jonas Ries, EMBL, jonas.ries@embl.de

% results_ax1=initaxis(par.resultstabgroup,'cross-correlations');
% results_ax3=initaxis(par.resultstabgroup,'normalized cross-correlations');

%here: rather from par, from channel range. Make sure it does not get
%displaced. Fill in outside.


lastframe=round(par.framestop);
firstframe=round(par.framestart);
numframes=lastframe-firstframe+1;

%% calculate movie and FFT of movie
pixrec=par.drift_pixrecz; %in nm
zb=par.zrange(1):pixrec:par.zrange(2);

window=ceil((par.drift_windowz-1)/2);
timepoints=par.drift_timepointsz; %how many timepoints
% maxdrift=par.drift_maxdrift; %in nanometers

mx=[min(pos.xnm) max(pos.xnm)]; %ROI which is used for drift correction. 
xb=mx(1):par.slicewidth:mx(2);
% my=[min(pos.ynm) max(pos.ynm)]; %You can put your own routine here

% srec(1)=round((mx(2)-mx(1))/pixrec);
% srec(2)=round((my(2)-my(1))/pixrec);

% srim= histrender(posr,mx, my, pixrec, pixrec);

% nfftexp=2^ceil(log2(max(srec))); %for fft use power of 2
% noff=nfftexp/2+1; 
% disp('make movie')
% Fmovier=makemovie;  %calculate fourier transforms of reconstructed images
% disp('find displacement')
plotaxis=initaxis(par.resultstabgroup,'CC z');
[ddz,errz]= finddisplacementsZ; %determine displacements
% 
% ddx=ddx*pixrec; %convert displacements into nm
% ddy=ddy*pixrec;

%% bin displacements
[dz,sdzc]=bindisplacementfit(ddz,errz); %determine displacement for each time point
% [dy,sdyc]=bindisplacementfit(ddy,erry);

% dx0h=[0; dx];
% dy0h=[0; dy];
s=size(ddz);
% ddxplot=cumsum(diff(ddx));
% ddyplot=cumsum(diff(ddy));
ddzplot=ddz;
% ddyplot=ddy;
for kn=1:s(1)
    ddzplot(:,kn)=ddz(:,kn)-ddz(kn,kn)+dz(kn);
%     ddyplot(:,kn)=ddy(:,kn)-ddy(kn,kn)+dy(kn);
end




%interpolate displacemnt for all frames
cfit1=(0:length(dz)-1)*binframes+binframes/2+firstframe; %positions of time points
ctrue=(1:par.maxframeall)'; %positions of frames
if 0%length(dz)>9
    [~,sdz,inlier,outlier]=robustMean(ddzplot,2,15);%std for each time point, used for interpolation
%     [~,sdy]=robustMean(ddyplot,2,15);
    sdzm=robustMean(sdz)/2;
    sdz(sdz<sdzm)=sdzm;

    sdzm=robustMean(sdz);

    indgz=sdz<10*sdzm;
%     if length(indgz)<9
%         indgz=true(size(dz));
%     end
else
    sdz=std(ddzplot,0,2); %std for each time point, used for interpolation

%     indgz=true(size(dz));
end

   indgz=true(size(dz));


% sdx(sdx>5*sdxm)=inf;
% sdy(sdy>5*sdym)=inf;


wz=1./sdz.^2;


%add to weights
% wx=wx+0.0*mean(wx);
% wy=wy+0.0*mean(wy);

% wx(end)=wx(end)/4;
% wy(end)=wy(end)/4;

 h=cfit1(2)-cfit1(1);
 if ~isempty(par.smoothpar)
    pset=1/(1+h^3/60*par.smoothpar);
 else
     pset=[];
 end
%  pset=1/(1+h^3/24);
%give higher weight to first data point:
% wx(1)=wx(1)*2;
% wy(1)=wy(1)*2;
% pset=[];

switch par.smoothmode.Value
    case 1 % smoothing spline
        [dzt,pz] = csaps(double(cfit1(indgz)),double(dz(indgz)),double(pset),double(ctrue),wz(indgz)) ;    
    case 2
        dzt = interp1(double(cfit1(indgz)),double(dz(indgz)),double(ctrue)) ;
end

% [dzt,pz] = csaps(double(cfit1(indgz)),double(dz(indgz)),pset,double(ctrue),wz(indgz)) ;

framesall=(1:par.maxframeall);%-firstframe+1;
binend=floor(1*binframes/2);
% dxtt=zeros((par.maxframeall),1);dytt=dxtt;
dztt=dzt;
dztt(1:firstframe-1+binframes/2)=dztt(firstframe-1+binframes/2+1);
% dxtt(firstframe:lastframe)=dxt;
dztt(lastframe+1-binend:end)=dztt(lastframe+1-binend);


results_ax2=initaxis(par.resultstabgroup,'dz/frame');

% subplot(1,2,1)
hold off
plot(ddzplot)
hold on
plot(dz,'k','LineWidth',1.5);
plot(sdz,'k:')
sz=(max(dz)-min(dz));
ylim([min(dz)-sz/2 max(dz)+sz/2])
axis tight


if par.drift_reference
    dztt=dztt-dz(end-1);

    dz=dz-dz(end-1);
end
    

driftinfo.dz=dz;

driftinfo.dzplot=ddzplot;

driftinfo.dzt=dztt;

driftinfo.binframesz=cfit1;
%
initaxis(par.resultstabgroup,'dz/frame final');

% hold off
plot(cfit1,dz,'x',framesall,dztt,'k')
% hold on
% plot(cfit1,dy,'o',framesall,dytt,'r')
xlabel('frame')
ylabel('dz (nm)')
drawnow



drift.z=dztt;


% asdafd
% fitposc=adddrift(positions,dxt,dyt); %recalculate positions



function [zpos,errz]= finddisplacementsZ % find displacements
binframes=2*ceil(numframes/timepoints/2+1);
frameranges=[firstframe:binframes:lastframe lastframe] ;   
dnumframesh =length(frameranges);
zpos=zeros(dnumframesh-1);
errz=zpos+1;

timerh=tic;
for k=1:dnumframesh-1
    indframek=pos.frame<frameranges(k+1)&pos.frame>=frameranges(k);
%     posr.x=pos.xnm(indframe);posr.y=pos.ynm(indframe);
    for l=k+1:dnumframesh-1
        if toc(timerh)>.5
            drawnow
            timerh=tic;
            plotaxish=plotaxis;
        else
            plotaxish=[];
        end
        indframel=pos.frame<frameranges(l+1)&pos.frame>=frameranges(l);
        zpos(k,l)=finddisplacementZ(pos.xnm(indframek),pos.znm(indframek),pos.xnm(indframel),pos.znm(indframel),xb,zb,window, plotaxish);
        zpos(l,k)=-zpos(k,l);
        
        if SMAP_stopnow
            error('execution stopped by user');
        end
   
    end
%     disp(k/dnumframesh)
end
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
