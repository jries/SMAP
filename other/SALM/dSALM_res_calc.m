function [zxdSALM_err]=dSALM_res_calc(zax,phot,bg, bgsfac, Ntotconst, calibfile)
%calculate dSALM error based on input

%global parameters

bgs=bg*bgsfac;
pixelsize=117;

roisize=59;

% Ntotconst=false; %if true, the total intensiyt (SAF+UAF) is constant, valid for fluorophore with high quantum yield. Otherwise UAF is constant. Valied for fluorophore with low quantum yield.
% userJe=false; %if to use SAF/UAF ratio from Jesacher PSF. Otherwise calculate with my funtion
% zax=-500:10:500;
% dz=zax(2)-zax(1);
% indplot=51:101;

coord=zeros(length(zax),3);coord(:,3)=-zax; %where to evaluate PSFs. 

%% intensities
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

%% our experimental PSF
% compare global calibration with normalization (fix)

psfglobfile=calibfile;
pixeliszeexp=117;

psfxuaf=splinePSF;psfxuaf.loadmodel(psfglobfile,1);
psfxsaf=splinePSF;psfxsaf.loadmodel(psfglobfile,2);

[CXsaf,CXsafP]=(psfxsaf.crlb(Ns_17,bgs,coord*0,roisize));
[CXuaf,CXuafP]=(psfxuaf.crlb(Nu_17,bgs,coord,roisize));
Nuafx_err=sqrt(CXuafP(:,1));
Nsafx_err=sqrt(CXsafP(:,1));

[zxdSALM_err]=zerr(Ns_17, Nu_17,Nsafx_err,Nuafx_err,dz);

% defocusexp= 200;
% CXsuafdef=(psfxuaf.crlb(Nu_17,bg,zax'+defocusexp,roisize));
% zxsdSALM_combined_err=1./sqrt(1./CXsuafdef(:,5)+1./zxdSALM_err.^2);
% xerrdexp=sqrt(CXsuafdef(:,1))*pixeliszeexp;
% yerrdexp=sqrt(CXsuafdef(:,2))*pixeliszeexp;
% avxerrdSALM=(xerrdexp+yerrdexp)/2;

% calculate vSALM based on uaf PSF fro NA1.7
% Nuv17=Nu_17/2;
% Ntotv17=(Nu_17+Ns_17)/2;
% [CXtot17,CXtot17P]=(psfxuaf.crlb(Ntotv17,(bgs+bg)/2,coord,roisize));
% [CXuaf17,CXuaf17P]=(psfxuaf.crlb(Nuv17,bg/2,coord,roisize));
% Nuaf17_err=sqrt(CXuaf17P(:,1));
% Ntot17_err=sqrt(CXtot17P(:,1));
% [zvSALM17_err]=zerr(Ntotv17, Nuv17,Ntot17_err,Nuaf17_err,dz);


end
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

