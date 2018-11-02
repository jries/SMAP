function [drift,driftinfo]=finddriftfeatureM(pos,par)
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

%other functions needed:
%myhist2
%my2DGaussfit

%copyright: Jonas Ries, EMBL, jonas.ries@embl.de

results_ax1=initaxis(par.resultstabgroup,'CC');
axiscurve=initaxis(par.resultstabgroup,'CC vs M');

%here: rather from par, from channel range. Make sure it does not get
%displaced. Fill in outside.


lastframe=round(par.framestop);
firstframe=round(par.framestart);
numframes=lastframe-firstframe+1;

%% calculate movie and FFT of movie
pixrec=par.drift_pixrec; %in nm
dM=((par.drift_window));
timepoints=par.drift_timepoints; %how many timepoints
% maxdrift=par.drift_maxdrift; %in nanometers

mx=[min(pos.xnm) max(pos.xnm)]; %ROI which is used for drift correction. 
my=[min(pos.ynm) max(pos.ynm)]; %You can put your own routine here

srec(1)=round((mx(2)-mx(1))/pixrec);
srec(2)=round((my(2)-my(1))/pixrec);

% srim= histrender(posr,mx, my, pixrec, pixrec);

nfftexp=2^ceil(log2(max(max(srec),256))); %for fft use power of 2
% noff=nfftexp/2+1; 
disp('make movie')

dummyim=histrender(struct('x',mean(mx),'y',mean(my)),mx, my, pixrec, pixrec)';
sdim=size(dummyim);
midp=round(sdim/2);
nx=((1:sdim(1))-midp(1))/sdim(1);
ny=((1:sdim(2))-midp(2))/sdim(2);
[X,Y]=meshgrid(nx,ny);
maskimage=sqrt(X.^2+Y.^2);


[Fmovier,posall]=makemovie;  %calculate fourier transforms of reconstructed images
disp('find displacement')
[Mall,Aall]= finddisplacements2; %determine displacements
[Mrel,sx]=bindisplacementfit(Mall,Aall);
mfit=cumprod([1; Mrel]);


Malln=Mall;
Malln(Malln==1)=NaN;
for k=1:size(Mall,1)
    Malln(k,:)=Malln(k,:)*mfit(k);
end
% Malln=triu(Malln);



Mfit=bindisp(Mrel);
Mfitn=Mfit;
Mfitn(Mfitn==1)=NaN;
for k=1:size(Mfit,1)
    Mfitn(k,:)=Mfitn(k,:)*mfit(k);
end





%interpolate displacemnt for all frames
cfit1=(0:length(mfit)-1)*binframes+binframes/2+firstframe; %positions of time points
ctrue=(1:par.maxframeall)'; %positions of frames



[mxt,pmfit] = csaps(double(cfit1),double(mfit),[],double(ctrue)) ;


% framesall=(1:par.maxframeall)-firstframe+1;
binend=floor(1*binframes/2);
binend=1;
% dxtt=zeros((par.maxframeall),1);dytt=dxtt;
% dxtt=dxt;
% mxt(1:firstframe-1+binframes/2)=dxtt(firstframe-1+binframes/2);
% % dxtt(firstframe:lastframe)=dxt;
mxt(lastframe+1-binend:end)=mxt(lastframe+1-binend);
% 
% dytt=dyt;
% dytt(1:firstframe-1+binframes/2)=dytt(firstframe-1+binframes/2);
% % dytt(firstframe:lastframe)=dyt;
% dytt(lastframe+1-binend:end)=dytt(lastframe+1-binend);


results_ax2=initaxis(par.resultstabgroup,'dM/frame');

plot(cfit1,Malln','-x');
hold on;
% plot(cfit1,mfit,'k');hold off;
% results_ax3=initaxis(par.resultstabgroup,'dM fit');


% hold on
plot(ctrue,mxt,'r','LineWidth',2)
plot(cfit1,mfit,'k*')
hold off
drawnow
if par.drift_reference
    mxt=mxt/mfit(end);
end
    

driftinfo.mframe=mxt;

driftinfo.mblock=mfit;

driftinfo.Mmat=Malln';
% driftinfo.dyt=dytt;
driftinfo.binframes=cfit1;

drift.M=mxt;



function [Fmovier,posh]=makemovie %calculate fourier transforms of images
    
%     posr.x=pos.xnm;posr.y=pos.ynm;
    binframes=2*ceil(numframes/timepoints/2+1);
    frameranges=[firstframe:binframes:lastframe lastframe] ;  
    Fmovier=zeros(nfftexp,nfftexp,timepoints,'single');
%     imager=zeros(nfftexp,nfftexp,timepoints,'single');
    for k=timepoints:-1:1
        indframe=pos.frame<frameranges(k+1)&pos.frame>=frameranges(k);
        posr.x=pos.xnm(indframe);posr.y=pos.ynm(indframe);
        imager(:,:,k)=(histrender(posr,mx, my, pixrec, pixrec).*maskimage)';
        Fmovier(:,:,k)=fft2(imager(:,:,k),nfftexp,nfftexp);
        posh(k)=posr;
        
        if SMAP_stopnow
            error('execution stopped by user');
        end
%         figure(89)
%         imagesc(imager)
%         waitforbuttonpress
    end
end

function [Mall,Aall]= finddisplacements2 % find displacements
% global SMAP_stopnow;
    %     dM=0.006;
s=size(Fmovier);
dnumframesh =s(3);
Mall=ones(dnumframesh-1);
Aall=eye(dnumframesh-1);
% fhold=imagesc(1,'Parent',results_ax1);

for k=1:dnumframesh-1
    oldmag=1;
    for l=k+1:dnumframesh
        oldmag=getMag(posall(l),Fmovier(:,:,k),oldmag,dM);
        title(axiscurve,[k; l])
        [oldmag,Acc]=getMag(posall(l),Fmovier(:,:,k),oldmag,dM/6);
        Mall(k,l)=oldmag;
        Aall(k,l)=Acc;
%         if isfield(par,'showresults') && par.showresults && toc(timerh)>0.5
%             timerh=tic;
%             fhold=imagesc(outim,'Parent',results_ax1);
%             imagesc(outimnorm,'Parent',results_ax3)
%             results_ax3.Title.String=num2str(k/dnumframesh+(l-k)/dnumframesh^2);
%             results_ax1.Title.String=num2str(k/dnumframesh+(l-k)/dnumframesh^2);
%             drawnow
%         end
    title(axiscurve,[k; l])
        if SMAP_stopnow
            error('execution stopped by user');
        end
    end
    

%     disp(k/dnumframesh)
end
end

function [mag,ampl]=getMag(pos,im2,magnhere,dM)

    win=4;
% dMa=[-2*dM;-dM;0;dM;2*dM];
% dMa=[-dM;0;dM];
dMa=[-3*dM;-2*dM;-dM;0;dM;2*dM;3*dM];
 [maxv]=getmaxv(pos,im2,magnhere,dMa);   
 magt=magnhere+dMa;
%   magt=(magnif+magnhere);
fp=fit(magt,maxv,'poly2');
newmag=-fp.p2/2/fp.p1;
magt2=magt(1):dM/15:magt(end);
plot(axiscurve,magt,maxv,'*',magt2,fp(magt2),newmag,fp(newmag),'ko');drawnow

% mag=magnhere+newmag-1;  
mag=newmag;
ampl=fp(newmag);
    function [maxv]=getmaxv(pos,im2,magnhere,dMa)
       
        
        maxv=zeros(3,1);
%          matr=zeros(3);
        for M=1:length(dMa)
            posh.x=(pos.x-mean(mx))*(dMa(M)+magnhere)+mean(mx);
            posh.y=(pos.y-mean(my))*(dMa(M)+magnhere)+mean(my);
            imager=(histrender(posh,mx, my, pixrec, pixrec).*maskimage)';
            fmovt=fft2(imager,nfftexp,nfftexp);
%            matr(1,1)=dMa(M)+magnhere;matr(2,2)=dMa(M)+magnhere;matr(3,3)=1;
%             tf=affine2d(matr);
%             R=imref2d(size(im1));
%             movt=imwarp(im1,tf,'OutputView',R);
%             fmovt=fft2(movt,nfftexp,nfftexp);
            cc=fmovt.*conj(im2);
            ccf=fftshift(ifft2(cc));
            s=round(size(ccf)/2);
            ssm=[100;100];
            intim=ccf(s(1)-ssm(1):s(1)+ssm(1),s(1)-ssm(1):s(1)+ssm(1));
            h=fspecial('gaussian',5,1.5);
             intimf=filter2(h,intim);
        
            [maxvi,indm]=(max(intimf(:)));
             
            [xm,ym]=ind2sub(size(intimf),indm);
            smallim=intimf(xm-win:xm+win,ym-win:ym+win);

            imbig=imresize(smallim,3,'bicubic');
            maxv(M)=max(imbig(:));%/(dMa(M)+magnhere);
           imagesc(results_ax1,intim);title(results_ax1,maxv(M));drawnow    
        end
        [~,maxind]=max(maxv);
        if maxind<=2||maxind>=length(maxv)-1 %out of range
            
%             figure(89);plot(magnhere+dMa,maxv)
            magnhere=magnhere+dMa(maxind);
            maxv=getmaxv(pos,im2,magnhere,dMa);

        end
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

function [fp,sx]=bindisplacementfit(M,A)

%idea: we measure displacements between every frames (dxik=xi-xk, xi is 
%displacement for frame i). Use all xi as fit parameters, fit function
%calculates dxik. Robust fit.
% weights=1./(errx+.1).^2;
%startp
sf=size(M);
% fp0=M(1,:);
fp0=ones(sf(1),1);
% fp0=zeros(sf(1)-1,1);
options=statset('nlinfit');
options=statset(options,'RobustWgtFun','');
A=A+eye(size(A))*max(A(:));
[fp,r,J,COVB,mse] = nlinfit(M(:),M(:),@bindispf,fp0,options,'Weights',A(:)+.1);
% options=statset(options,'RobustWgtFun','bisquare');
% [fp,r,J,COVB,mse] = nlinfit(M(:),M(:),@bindispf,fp,options,'Weights',A(:)+.1);
ci = nlparci(fp,r,'covar',COVB);
%  dx2=[0; fp];
%  dx2=fp;
sx=ci(:,2)-ci(:,1);
% sdx2=[mean(sx);sx];
end

function out=bindisp(fp,M)
lfp=length(fp);
out=ones(lfp,lfp+1);
for k=1:lfp
    out(k,k+1:end)=cumprod(fp(k:end));
end

% fph=[0; fp];
% ddxf=zeros(length(fph));
% for k=1:length(fph)
%     ddxf(k,:)=fph(k)-fph; %calculate difference
% end
% 
% 
% % out=ddxf(:);
% out=ddxf;
end
function out=bindispf(fp,M)
out=bindisp(fp,M);
out=out(:);
end
