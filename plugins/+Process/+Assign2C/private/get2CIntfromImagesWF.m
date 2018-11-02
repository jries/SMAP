function locout=get2CIntfromImagesWF(loc,transform,filestruc,p)
pix_cam=filestruc.info.cam_pixelsize_um*1000;
x=double(loc.xnm);
y=double(loc.ynm);

%clear a, b separation from beginning. a=ref side, b=target side

indref=transform.getRef(x,y);
indtarget=~indref;

xA=zeros(size(x));xB=zeros(size(x));yA=zeros(size(x));yB=zeros(size(x));

xA(indref)=x(indref);yA(indref)=y(indref);
[xA(indtarget),yA(indtarget)]=transform.transformCoordinatesInv(x(indtarget),y(indtarget));

xB(indtarget)=x(indtarget);yB(indtarget)=y(indtarget);
[xB(indref),yB(indref)]=transform.transformCoordinatesFwd(x(indref),y(indref));

fieldsout={};
%initialize
if p.checksum
    fieldsout={fieldsout{:},'sum_bg1','sum_bg2','sum_n1','sum_n2'};
end
if p.checkfit
    fieldsout={fieldsout{:},'fit_bg1','fit_bg2','fit_n1','fit_n2'};
end

for k=1:length(fieldsout)
    locout.(fieldsout{k})=zeros(length(x),1,'single');
end
tau=p.filtert;

path=fileparts(filestruc.name); %rather in top class, pass on
[~, path]=uigetfile([path filesep '*.tif'],filestruc.name);
df=myfastdir(path,  'img*.tif');
pfile.filelist=df;
pfile.path=path;
pfile.bufferSize=2*tau+1;



fileloader=fitter.FileLoader;
fileloader.setParameters(pfile);
pcut.kernelSize=13;
roiCutter=fitter.RoiCutter;
roiCutter.setParameters(pcut);

pbg.bgon=1;
pbg.rho=p.filterx;
pbg.cam_offset=filestruc.info.offset;
pbg.pix2phot=filestruc.info.pix2phot;

makeBG=fitter.MedianBackground;
makeBG.setParameters(pbg);


lf=length(loc.frame);
bgf=1;
ind2=1;


%determine PSFx and PSFy
if isfield(loc,'znm')
    if isfield(filestruc.info,'fit') %determine from fitz
        fitz=filestruc.info.fit.fitzParameters;
        parx= [fitz(7) fitz(1) fitz(2) fitz(4) fitz(6) 0];
        pary= [fitz(7) fitz(8) fitz(3) fitz(5) -fitz(6) 0];
        PSFxpix=sigmafromz(pary,loc.znm/1000,1);
        PSFypix=sigmafromz(parx,loc.znm/1000,1);
    else
        d=0.42;
        g=-0.2;
        sx0=1.1;
        PSFxpix=sigmafromz_simple(loc.znm/1000,[d -g sx0]);
        PSFypix=sigmafromz_simple(loc.znm/1000,[d g sx0]);
    end
else
    if p.fitUsePSF
        PSFxpix=double(loc.PSFxnm/pix_cam(1));
    else
        PSFxpix=ones(length(loc.xnm),1)*p.PSFxnm/pix_cam(1);
    end 
    PSFypix=PSFxpix;
end


for frame=min(loc.frame):max(loc.frame)
    if bgf<=frame %make new background
        range=frame:frame+tau-1;
        bg=adu2phot(makeBG.run(fileloader.getBlock(range)),pbg);
        bgmean=mean(bg(~isnan(bg)));
        bgf=frame+tau;
        disp(['frame: ' num2str(frame)])
    end
    
    
    
% compare intensities
    %find indices for same frame
    ind1=ind2;
    while loc.frame(ind1)<frame && ind1<lf;
        ind1=ind1+1;
    end
    ind2=ind1;
    while loc.frame(ind2)==frame && ind2<lf;
        ind2=ind2+1;
    end
    imageframe=adu2phot(fileloader.get(frame),pbg);
    
    [maxposA,maxposAR]=nm2pixLoc(xA(ind1:ind2),yA(ind1:ind2),pix_cam,filestruc.info.roi);
    roiA=roiCutter.run(imageframe,maxposAR);
    bgA=roiCutter.run(bg,maxposAR);
    
    [maxposB,maxposBR]=nm2pixLoc(xB(ind1:ind2),yB(ind1:ind2),pix_cam,filestruc.info.roi);
    roiB=roiCutter.run(imageframe,maxposBR);
    bgB=roiCutter.run(bg,maxposBR);
    
%   sum;
    locsum=roi2int_sum(roiA, roiB,bgA,bgB,p.sizeRoiSum);
%     locsum=assignside(locroi,{'sum_n','sum_bg'},x(ind1:ind2),y(ind1:ind2),transform.refPart,midp);
    

    
    %fit
    if ~p.fitOnBg
        bgAf=[];
        bgBf=[];
    else
        bgAf=bgA;
        bgBf=bgB;
    end
            %fit on BG, don't fit bg
        pA=roi2int_fit_e(roiA,maxposA.x-maxposAR.x,maxposA.y-maxposAR.y,PSFxpix(ind1:ind2),PSFypix(ind1:ind2),p.sizeRoiFit,bgAf);
        locfit.fit_n1=pA(:,1);
        locfit.fit_bg1=pA(:,2);  
        pB=roi2int_fit_e(roiB,maxposB.x-maxposBR.x,maxposB.y-maxposBR.y,PSFxpix(ind1:ind2),PSFypix(ind1:ind2),p.sizeRoiFit,bgBf);
        locfit.fit_n2=pB(:,1);
        locfit.fit_bg2=pB(:,2);
%     else
%         
%         pref=roi2int_fit_e(roiA,maxposA.x-maxposAR.x,maxposA.y-maxposAR.y,PSFxpix(ind1:ind2),PSFypix(ind1:ind2),p.sizeRoiFit);
%         locfit.fit_n1=pref(:,1);
%         locfit.fit_bg1=pref(:,2);  
%         ptarget=roi2int_fit_e(roitarget,maxposB.x-maxposBR.x,maxposB.y-maxposBR.y,PSFxpix(ind1:ind2),PSFypix(ind1:ind2),p.sizeRoiFit);
%         locfit.fit_n2=ptarget(:,1);
%         locfit.fit_bg2=ptarget(:,2);
%     end
%     locfitsort=assignside(locfit,{'fit_n','fit_bg'},x(ind1:ind2),y(ind1:ind2),transform.refPart,midp);
    
    
    locall=copyfields(locsum,locfit);

    for k=1:length(fieldsout)
        locout.(fieldsout{k})(ind1:ind2)=locall.(fieldsout{k});
    end
end

function PSFx=sigmafromz_simple(z,p)%[d g sx0]);
    PSFx=p(3).*sqrt(1+(z-p(2)).^2./p(1).^2);



function s=sigmafromz(par,z,B0)
par=real(par);
% parx= [d sx0 Ax Bx g mp]
s0=par(2);d=par(1);A=par(3);B=par(4)*B0;g=par(5);mp=par(6);

% s=s0*sqrt(1+(z-g+mp).^2/d^2);
s=s0*sqrt(1+(z-g+mp).^2/d^2+A*(z-g+mp).^3/d^3+B*(z-g+mp).^4/d^4);
s=real(s);


%     %                    setbusy(par.handles.status,['gpu fit. frame: ' num2str(fr-par.firstframe) ' of ' num2str(length(frames)) ...
% %                        ' (' num2str((fr-par.firstframe+1)/(par.lastframe-par.firstframe+1)*100,'%10.0f')...
% %                        '%), , time = ' num2str(tv,'%3.0f') ' of ' num2str(tall,'%3.0f') 's, locs=' num2str(loctot,'%3.0f')])
% %                    drawnow


% figure(88);
% plot(loc.xnm,loc.ynm,'r.');
% hold on
% plot(xo,yo,'b.');
% hold off
% p.dataselect.value




function out=assignside(in,fieldpref,x,y,direction,midp)
switch direction
    case 'top-bottom'
        indleft=y<midp;
    case 'left-right'
        indleft=x<midp;
    otherwise
        disp('all transform not implemented')
end
    indright=~indleft;

for k=1:length(fieldpref)
    out.([fieldpref{k} '1'])(indleft)=in.([fieldpref{k} '1'])(indleft);
    out.([fieldpref{k} '1'])(indright)=in.([fieldpref{k} '2'])(indright);
    out.([fieldpref{k} '2'])(indleft)=in.([fieldpref{k} '2'])(indleft);
    out.([fieldpref{k} '2'])(indright)=in.([fieldpref{k} '1'])(indright);    
end

% function [loc,locr]=nm2pixLoc(x,y,pixelsize,roi)
% loc.x=(x/pixelsize(1))-roi(1);
% loc.y=(y/pixelsize(2))-roi(2);
% locr.x=round(loc.x);
% locr.y=round(loc.y);

function out=adu2phot(in,p)
out=(in-p.cam_offset)*p.pix2phot;




% function [xo,yo]=applytform(x,y,transform,midp)
% x=double(x);y=double(y);
% switch transform.refpart
%     case 2 %left
%         ind1=x<=midp;
%     case 3%right
%         ind1=x>=midp;
%     case 4 %top
%         ind1=y<=midp;
%     case 5% bottom
%         ind1=y>=midp;
%     otherwise %all
%         dips('only works for two parts on chip')
% end
% xo=zeros(size(x));
% yo=zeros(size(y));
% [xo(ind1),yo(ind1)]=transformPointsInverse(transform.T,x(ind1),y(ind1));
% [xo(~ind1),yo(~ind1)]=transformPointsForward(transform.T,x(~ind1),y(~ind1));