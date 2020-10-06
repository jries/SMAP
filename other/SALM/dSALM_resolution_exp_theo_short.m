%based on compareCRBSplineJesacher
%checked, annotated, simplified for dSALM paper

%global parameters
phot=5000;
bg=50;
bgs=bg/10;
pixelsize=117;
roisize=59;
Ntotconst=false; %if true, the total intensiyt (SAF+UAF) is constant, valid for fluorophore with high quantum yield. Otherwise UAF is constant. Valied for fluorophore with low quantum yield.
userJe=false; %if to use SAF/UAF ratio from Jesacher PSF. Otherwise calculate with my funtion
zax=-500:10:500;
dz=zax(2)-zax(1);
indplot=51:101;

coord=zeros(length(zax),3);coord(:,3)=-zax; %where to evaluate PSFs. 


%% Theoretical PSFs
%these were caluclated with the software that was provided with the
%Jesacher paper
crlpath='/Users/jonasries/Documents/MATLAB/other/CRB_calc/';
addpath([crlpath 'classes'],[crlpath 'functions']);
psfpath='example_PSFs/';

%Olympus 1.7
PSFuaf_tJ17='PSF_def0_(uaf)_NA1.7_72x72_-500-10-500nm_RI=1.33';
PSFdef_tJ17='PSF_def-5e-07_(tot)_NA1.7_72x72_-500-10-500nm_RI=1.33';
PSFsaf_tJ17=strrep(PSFuaf_tJ17,'_(uaf)_','_(saf)_');
PSFtot_tJ17=strrep(PSFuaf_tJ17,'_(uaf)_','_(tot)_');

%the Olympus 100x 1.49
PSFuaf_tJ49='PSF_def0_(uaf)_NA1.49_72x72_-500-10-500nm_RI=1.33';
PSFtot_tJ49=strrep(PSFuaf_tJ49,'_(uaf)_','_(tot)_');
PSFsaf_tJ49='PSF2_d0_(saf)_NA1.49_72x72_-500-10-500nm_RI=1.33';
PSFsaf_tJ49='PSF_os1_d0_(saf)_NA1.49_72x72_-500-10-500nm_RI=1.33';
%load PSF models
% spline PSFs
[psfSuaf_tJS17,normuaf_tJS17l]=loadspline([crlpath psfpath],PSFuaf_tJ17,roisize);
[psfSsaf_tJS17,normsaf_tJS17l]=loadspline([crlpath psfpath],PSFsaf_tJ17,roisize);
[psfStot_tJS17,normtot_tJS17l]=loadspline([crlpath psfpath],PSFtot_tJ17,roisize);

[psfSuaf_tJS49,normuaf_tJS49l]=loadspline([crlpath psfpath],PSFuaf_tJ49,roisize);
[psfStot_tJS49,normtot_tJS49l]=loadspline([crlpath psfpath],PSFtot_tJ49,roisize);
[psfSsaf_tJS49,normsaf_tJS49l]=loadspline([crlpath psfpath],PSFsaf_tJ49,roisize);
% 
% %calculate image stacks
% psfimSuaf_tJ17=double(psfSuaf_tJS17.PSF(coord));
% psfimSsaf_tJ17=double(psfSsaf_tJS17.PSF(coord));
% psfimStot_tJ17=double(psfStot_tJS17.PSF(coord));

%normalized to 1, now follow Je and normalize to Nsaf+Nuaf=Ntot
normsaf_tJS17=normsaf_tJS17l./(normsaf_tJS17l+normuaf_tJS17l).*normtot_tJS17l;
normuaf_tJS17=normuaf_tJS17l./(normsaf_tJS17l+normuaf_tJS17l).*normtot_tJS17l;
normtot_tJS17=normtot_tJS17l;

normsaf_tJS49=normsaf_tJS49l./(normsaf_tJS17l+normuaf_tJS17l).*normtot_tJS17l;
normuaf_tJS49=normuaf_tJS49l./(normsaf_tJS17l+normuaf_tJS17l).*normtot_tJS17l;
normtot_tJS49=normtot_tJS49l;
%% CRLBs
%Spline
photStot_tJS17=phot.*normtot_tJS17/psfStot_tJS17.normalization;photSuaf_tJS17=phot.*normuaf_tJS17/psfSuaf_tJS17.normalization;photSsaf_tJS17=phot.*normsaf_tJS17/psfSsaf_tJS17.normalization;
[CStot_tJS17,CStotP_tJS17]=psfStot_tJS17.crlb(photStot_tJS17,bg,coord);
[CSuaf_tJS17,CSuafP_tJS17]=psfSuaf_tJS17.crlb(photSuaf_tJS17,bg,coord);[CSsaf_tJS17,CSsafP_tJS17]=psfSsaf_tJS17.crlb(photSsaf_tJS17,bgs,coord);

%% photometry-based resolution dSALM
% calculate the relative intensity in teh SAF channel
pdSALM_17.n1=1.33;
pdSALM_17.n2=1.78;
pdSALM_17.lambda=670;%psfJtot_tJ17.lambda*1e9; %nm as in make PSF
pdSALM_17.NA=1.70;
pdSALM_17.NAmask=1.33;
[r_17,Is_17,Iu_17]=intensitySALM(zax,pdSALM_17); 

if Ntotconst
    Ns_17=r_17./(1+r_17)*phot;
    Nu_17=1./(1+r_17)*phot;
else
    Nu_17=phot+0*r_17;
    Ns_17=r_17*phot;
end
Ntot_17=Nu_17+Ns_17;

Nuafd_err_tJS17=sqrt(CSuafP_tJS17(:,1)); %error in photons. 
Nsafd_err_tJS17=sqrt(CSsafP_tJS17(:,1));
[zdSALM_err_tJS17,dz_drd_tJS17,dRd_tJS17]=zerr(Ns_17, Nu_17,Nsafd_err_tJS17,Nuafd_err_tJS17,dz);

%% photometry-based resolution vSALM
%recalculate with half the photon count (splitting)

Bgtot2=(bgs+bg)/2;
Bguaf2=bg/2;%splitting
Ntot2=(Ntot_17)/2; %splitting
Nuaf2=Nu_17/2;

[CStot2_tJS17v,CStotP2_tJS17v]=psfStot_tJS17.crlb(Ntot2,Bgtot2,coord);
[CSuaf2_tJS17v,CSuafP2_tJS17v]=psfSuaf_tJS17.crlb(Nuaf2,Bguaf2,coord);
Nuafv_err_tJS17v=sqrt(CSuafP2_tJS17v(:,1));
Ntotv_err_tJS17v=sqrt(CStotP2_tJS17v(:,1));
[zvSALM_err_tJS17v,dz_drv_tJS17v,dRv_tJS17v]=zerr(Ntot2, Nuaf2,Ntotv_err_tJS17v,Nuafv_err_tJS17v,dz);


%% combined resolution
zdcomb_tJS17=1./sqrt(1./zdSALM_err_tJS17.^2+1./CSuaf_tJS17(:,5));
zvcomb_tJS17v=1./sqrt(1./zvSALM_err_tJS17v.^2+1./CSuaf2_tJS17v(:,5)+1./CStot2_tJS17v(:,5));

%% x, y resolution 
errSsavX_tJS17=sqrt(CSsaf_tJS17(:,1))*pixelsize;
errSuavX_tJS17=sqrt(CSuaf_tJS17(:,1))*pixelsize;
errStotX_tJS17=sqrt(CStot_tJS17(:,1))*pixelsize;
errSuavX2_tJS17v=sqrt(CSuaf2_tJS17v(:,1))*pixelsize;
errStotX2_tJS17v=sqrt(CStot2_tJS17v(:,1))*pixelsize;
xdcomb_tJS17=1./sqrt(1./errSuavX_tJS17.^2+1./errSsavX_tJS17.^2);
xvcomb_tJS17v=1./sqrt(1./errSuavX2_tJS17v.^2+1./errStotX2_tJS17v.^2);


%% vSALM 1.49
p149.n1=1.33;
p149.n2=1.518;
p149.lambda=670; %nm as in make PSF
p149.NA=1.49;
p149.NAmask=1.33;
% Ntotconst=true;
[r49]=intensitySALM(zax,p149);  %compare
if Ntotconst
    Ns49=r49./(1+r49)*phot;
    Nu49=1./(1+r49)*phot;
else
    Nu49=phot;
    Ns49=r49*phot;
end
Ntot49=Nu49+Ns49;
Ntot249=Ntot49/2;
Nuaf249=Nu49/2;
Bgtot249=(bgs+bg)/2;
Bguaf249=bg/2;

[CStot249,CStotP249]=psfStot_tJS49.crlb(Ntot249,Bgtot249,coord);
[CSuaf249,CSuafP249]=psfSuaf_tJS49.crlb(Nuaf249,Bguaf249,coord);
Nuafv_err49=sqrt(CSuafP249(:,1));
Ntotv_err49=sqrt(CStotP249(:,1));
[zv49SALM_err]=zerr(Ntot249, Nuaf249,Ntotv_err49,Nuafv_err49,dz);
errSuavX249=sqrt(CSuaf249(:,1))*pixelsize;
errStotX249=sqrt(CStot249(:,1))*pixelsize;

zvcomb49=1./sqrt(1./zv49SALM_err.^2+1./CSuaf249(:,5)+1./CStot249(:,5));
xvcomb49=1./sqrt(1./errSuavX249.^2+1./errStotX249.^2);

% dSALM 1.49

photStot_tJS49=phot.*normtot_tJS49/psfStot_tJS49.normalization;photSuaf_tJS49=phot.*normuaf_tJS49/psfSuaf_tJS49.normalization;photSsaf_tJS49=phot.*normsaf_tJS49/psfSsaf_tJS49.normalization;
[CStot_tJS49,CStotP_tJS49]=psfStot_tJS49.crlb(photStot_tJS49,bg,coord);
[CSuaf_tJS49,CSuafP_tJS49]=psfSuaf_tJS49.crlb(photSuaf_tJS49,bg,coord);
[CSsaf_tJS49,CSsafP_tJS49]=psfSsaf_tJS49.crlb(photSsaf_tJS49,bgs,coord);

Nuafd_err_tJS49=sqrt(CSuafP_tJS49(:,1)); %error in photons. 
Nsafd_err_tJS49=sqrt(CSsafP_tJS49(:,1));
[zdSALM_err_tJS49,dz_drd_tJS49,dRd_tJS49]=zerr(Ns49, Nu49,Nsafd_err_tJS49,Nuafd_err_tJS49,dz);


%% our experimental PSF
% compare global calibration with normalization (fix)
psfglobfile='/Volumes/LaCie/otherData/SALM/CRLB/dSALM/dSALM_global_3dcal.mat';
psfsaffile='/Volumes/LaCie/otherData/SALM/CRLB/dSALM/dSALM_saf_3dcal.mat';
psfuaffile='/Volumes/LaCie/otherData/SALM/CRLB/dSALM/dSALM_uaf_3dcal.mat';
pixeliszeexp=117;

psfxuaf=splinePSF;psfxuaf.loadmodel(psfglobfile,1);
psfxsaf=splinePSF;psfxsaf.loadmodel(psfglobfile,2);

[CXsaf,CXsafP]=(psfxsaf.crlb(Ns_17,bgs,coord*0,roisize));
[CXuaf,CXuafP]=(psfxuaf.crlb(Nu_17,bgs,coord,roisize));
Nuafx_err=sqrt(CXuafP(:,1));
Nsafx_err=sqrt(CXsafP(:,1));

[zxdSALM_err]=zerr(Ns_17, Nu_17,Nsafx_err,Nuafx_err,dz);

defocusexp= 200;
CXsuafdef=(psfxuaf.crlb(Nu_17,bg,zax'+defocusexp,roisize));
zxsdSALM_combined_err=1./sqrt(1./CXsuafdef(:,5)+1./zxdSALM_err.^2);
xerrdexp=sqrt(CXsuafdef(:,1))*pixeliszeexp;
yerrdexp=sqrt(CXsuafdef(:,2))*pixeliszeexp;
avxerrdSALM=(xerrdexp+yerrdexp)/2;

% calculate vSALM based on uaf PSF fro NA1.7
Nuv17=Nu_17/2;
Ntotv17=(Nu_17+Ns_17)/2;
[CXtot17,CXtot17P]=(psfxuaf.crlb(Ntotv17,(bgs+bg)/2,coord,roisize));
[CXuaf17,CXuaf17P]=(psfxuaf.crlb(Nuv17,bg/2,coord,roisize));
Nuaf17_err=sqrt(CXuaf17P(:,1));
Ntot17_err=sqrt(CXtot17P(:,1));
[zvSALM17_err]=zerr(Ntotv17, Nuv17,Ntot17_err,Nuaf17_err,dz);


%% also vSALM
psfglobfile149='/Volumes/LaCie/otherData/SALM/CRLB/vSALM/vSALM_stack_3dcal.mat';
psfxuaf49=splinePSF;psfxuaf49.loadmodel(psfglobfile149,1);
psfxtot49=splinePSF;psfxtot49.loadmodel(psfglobfile149,2);

[CXtot249,CXtotP249]=psfxtot49.crlb(Ntot249,Bgtot249,coord);
[CXuaf249,CXuafP249]=psfxuaf49.crlb(Nuaf249,Bguaf249,coord);
Nxuafv_err49=sqrt(CXuafP249(:,1));
Nxtotv_err49=sqrt(CXtotP249(:,1));
[zxv49SALM_err]=zerr(Ntot249, Nuaf249,Nxtotv_err49,Nxuafv_err49,dz);

defocusexp49= -200;
CXtotdef49=psfxtot49.crlb(Ntot249,Bgtot249,zax'+defocusexp49,roisize);
CXuafdef49=psfxuaf49.crlb(Ntot249,Bgtot249,zax'+defocusexp49,roisize);

errXuavX249=sqrt(CXtotdef49(:,1))*pixelsize;
errXtotX249=sqrt(CXuafdef49(:,1))*pixelsize;

errXuavY249=sqrt(CXtotdef49(:,2))*pixelsize;
errXtotY249=sqrt(CXuafdef49(:,2))*pixelsize;

zvpsf49=1./sqrt(1./CXtotdef49(:,5)+1./CXuafdef49(:,5));
zvcomb49=1./sqrt(1./zxv49SALM_err.^2+1./zvpsf49.^2);
xvcomb49=1./sqrt(1./errXuavX249.^2+1./errXtotX249.^2);
yvcomb49=1./sqrt(1./errXuavY249.^2+1./errXtotY249.^2);

%% plots
%z,x experimental (mian figure): v 1.49, v1.7, d1.7

fx=figure(875);
fx.Name='experimental PSF';
hold off
plot(zax(indplot), zxdSALM_err(indplot),'r')
hold on
plot(zax(indplot), sqrt(CXsuafdef(indplot,5)),'k:')
plot(zax(indplot), zxsdSALM_combined_err(indplot),'r--')
plot(zax(indplot), zxv49SALM_err(indplot),'b')
% plot(zax(indplot), avxerrdSALM(indplot),'m')
plot(zax(indplot), xerrdexp(indplot),'k')
plot(zax(indplot), yerrdexp(indplot),'k')
plot(zax(indplot), zvSALM17_err(indplot),'b-.')
plot(zax(indplot),zdSALM_err_tJS49(indplot),'b--')

legend('dSALM 1.7','Astig','dSALM+AS','vSALM 1.49','dSALM x','dSALM y','vSALM 1.7','dSALM 1.49*','location','northwest')
xlabel('z (nm)')
ylabel('localization precision (nm)')
title(['experimental PSF. Ntotconst: ' num2str(Ntotconst)])
ylim([0 30])
xlim([0 400])

%z theoretical based on spline PSF and Jesacher PSF model Z
ft=figure(876);
ft.Name='theoretical PSF';
hold off
plot(zax(indplot),zdSALM_err_tJS17(indplot),'r')
hold on
plot(zax(indplot),zdSALM_err_tJS49(indplot),'b--')
plot(zax(indplot), sqrt(CSuaf_tJS17(indplot,5)),'k-.')
plot(zax(indplot),zdcomb_tJS17(indplot),'r--')
plot(zax(indplot),zv49SALM_err(indplot),'b')
plot(zax(indplot), errSuavX_tJS17(indplot),'k')
plot(zax(indplot), zvSALM_err_tJS17v(indplot),'b-.')
% plot(zax(indplot), sqrt(Cz_def(indplot)),'g--')


legend('dSALM 1.7','dSALM 1.49','PSF','dSALM+PSF','vSALM 1.49','dSALM x,y','vSALM 1.7','location','northwest')
xlabel('z (nm)')
ylabel('localization precision (nm)')
title(['calculated PSF (Jesacher). Ntotconst: ' num2str(Ntotconst)])
ylim([0 30])
xlim([0 400])

fr=figure(877);
fr.Name='improvement z err';

hold off
plot(zax(indplot),zvSALM_err_tJS17v(indplot)./zdSALM_err_tJS17(indplot),'k')
hold on
plot(zax(indplot),zv49SALM_err(indplot)./zdSALM_err_tJS17(indplot),'k--')

xlabel('z (nm)')
ylabel('ratio z err vSALM vs dSALM')
title('ratio z err vSALM vs dSALM')
legend('vSALM NA1.7','vSALM NA1.49')

figure(878)
zaxh=(0:pdSALM_17.lambda);
[r_17p,Is_17p,Iu_17p]=intensitySALM(zaxh,pdSALM_17); 
rplot=Is_17p./(Iu_17p+Is_17p);
plot(zaxh/pdSALM_17.lambda, rplot)

%%
function [psf,norm]=loadspline(path,PSFname,roisize)
    psf=splinePSF;
    psf.loadmodel([path PSFname '_3dcal.mat'])
    psf.roisize=roisize;
    normspline= load([path PSFname '_normalization.mat']);
    norm=squeeze(normspline.norm);
    norm=norm/max(norm(round(end/2):end)); %XXXX only positive z
end


function [zdSALM_err,dz_dr,dR]=zerr(Ntot, Nuaf,Ntot_err,Nuaf_err,dz)
dru=-Ntot./Nuaf.^2;
drt=1./Nuaf;
dR2=drt.^2.*Ntot_err.^2+dru.^2.*Nuaf_err.^2;
r=Ntot./Nuaf;
dz_dr=1./(diff(r)/dz);dz_dr(end+1)=dz_dr(end);
dz2=dz_dr.^2.*dR2;
dR=sqrt(dR2);
zdSALM_err=sqrt(dz2);
end

