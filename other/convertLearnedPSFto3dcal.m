% importLearnedPSF(filename)
%% load
% persistent pfad
if ~exist('pfad','var')
    pfad='';
end
[f,pfad]=uigetfile([pfad '*.h5']);
file=[pfad f];
%% import
p.dz=double(h5read(file,'/dz'));
PSF=h5read(file,'/learned_psf');
info=h5info(file);
if length(info.Datasets)>2
    rois=single(h5read(file,'/rois'));
else 
    rois=[];
end
% PSF=permute(PSF,[2 1 3]);
%% smoothing bspline
p.smoothz=1;
lambdaz=p.smoothz/p.dz*100;
lambda=[0 0 lambdaz];
b3_0r=bsarray(double(PSF),'lambda',lambda);
%% make cspline
zhd=1:1:b3_0r.dataSize(3);
dxxhd=1;
[XX,YY,ZZ]=meshgrid(1:dxxhd:b3_0r.dataSize(1),1:dxxhd:b3_0r.dataSize(2),zhd);
% p.status.String='calculating cspline coefficients in progress';drawnow
corrPSFhdr = interp3_0(b3_0r,XX,YY,ZZ,0);
%calculate cspline coefficients
coeffr = single(Spline3D_interp(corrPSFhdr));

%% test
midp=size(PSF,3)/2;
 rangef=1:141;
%  rangef=1:2*midp;
roisize=15;
dx=floor(roisize/2);
roiP=ceil(size(PSF,1)/2);
rangex=roiP-dx:roiP+dx;
PSFfit=single(PSF(rangex,rangex,rangef));
 [P,CRLB, LL] =mleFit_LM(PSFfit*1000+30,5,50,coeffr,0,0,0);
 figure(22);
 subplot(2,2,1)

 plot(rangef,P(:,1)-dx+1);hold on; 
 plot(rangef,P(:,2)-dx+1);hold off;
 
 xlabel('frame')
 ylabel('x, y (pixel)')
 subplot(2,2,2)
 plot(rangef,P(:,5)-rangef'+1);
  xlabel('frame')
 ylabel('zfit-frame (frame)')
 
 subplot(2,2,3)
 hold off 
 plot(rangef,rangef,'k:')
hold on

 subplot(2,2,4)
 hold off 
 plot(rangef,rangef*0,'k:')
hold on
 if ~isempty(rois)
     numb=size(rois,4);
     for k=1:numb
         roih=rois(rangex,rangex,rangef,k);
         [Pb,CRLB, LL] =mleFit_LM(roih-min(roih(:))+10,5,50,coeffr,0,0,0);
          subplot(2,2,3)
          plot(rangef,Pb(:,5));
           subplot(2,2,4)
          plot(rangef,Pb(:,5)-rangef'+1);
%           plot(P(:,6));
     end
 end
           subplot(2,2,3)
  plot(rangef,P(:,5),'k')
  xlabel('frame')
 ylabel('z (frame)')
 
            subplot(2,2,4)
  plot(rangef,P(:,5)-rangef'+1,'k')
  xlabel('frame')
 ylabel('zfit-frame (frame)')
 
%% assemble output structure
S.cspline.dz=single(dz);
S.cspline.z0=single(round(midp));
S.cspline.coeff=coeffr;
S.Xrange=[-inf inf];
S.Yrange=[-inf inf];
S.EMon=false;
SXY=S;

%% save _3dcal.mat
outf=strrep(file,'.h5','_3dcal.mat');
save(outf,'SXY')
