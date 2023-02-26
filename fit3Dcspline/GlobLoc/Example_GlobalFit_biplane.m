% Copyright (c) 2021 Li Lab, Southern University of Science and Technology, Shenzhen & Ries Lab, European Molecular Biology Laboratory, Heidelberg.
% author: Yiming Li
% email: liym2019@sustech.edu.cn
% date: 2021.08.27
% Download the calbration data from the followling link first: https://www.embl.de/download/ries/globLoc/
% BP-combine_3dcal.mat
% Tested with CUDA 11.3 (Express installation) and Matlab 2019a
%% Load and make cspline PSF
clearvars
cal=load('BP-combine_3dcal.mat'); % PSF model from SMLM2016 challenge

PSF1 = cal.SXY(1).PSF{1};
PSFnorm1 = PSF1/mean(sum(sum(PSF1(:,:,49:51),1,'omitnan'),2,'omitnan'),'omitnan'); % normalize to focal spots
coeff1 = single(Spline3D_interp(PSFnorm1));
dz1=10; % distance between slides (nm)
z01=76; % central slides

PSF2 = cal.SXY(2).PSF{1};
PSFnorm2 = PSF2/mean(sum(sum(PSF2(:,:,99:101),1,'omitnan'),2,'omitnan'),'omitnan'); % normalize to focal spots
coeff2 = single(Spline3D_interp(PSFnorm2));
dz2=10; % distance between slides (nm)
z02=76; % central slides

%% Parameters

Nfits = 1000;
Nphotons1 =2500;
Nphotons2 =2500;
Npixels = 15;
bg1 = 20;
bg2 = 20;
zoffset = 0;%nm

ztruthAll = -600:50:600;
numz=length(ztruthAll);

zstartnm = 200;
zstart=repmat([-2 -1 1 2]*zstartnm/dz1,[Nfits 1]);
theta = 1*pi/180;%rotation
sx = 1.01;%scalex
sy = 1.02;%scaley
tx = -.5;%shfitx
ty = -0.2;%shfity

%define affine transformation matrix
tformR = [cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0;0 0 1];
tformS = [sx 0 0;0 sy 0;0 0 1];
tformT = [1 0 tx; 0 1 ty; 0 0 1];
tformF = tformR*tformS*tformT;
tforminv=inv(tformF);

%% simulate and fit
fgood_linkALL=zeros(numz,1);
fgoodl_linkXYZ=zeros(numz,1);
fgood_sep=zeros(numz,1);
for iz =1:numz
    disp(iz)
    coordsxy1 = Npixels/2 -1 +2*rand([Nfits 2]); 
    coordsxy2 = zeros(Nfits,2);
    
    for i = 1:Nfits
        temp = tformF*[coordsxy1(i,:)';1];
        coordsxy2(i,:) = temp(1:2); 
    end
    ztruth = ztruthAll(iz);
    dzrand=(rand(Nfits,1)-0.5)*(ztruthAll(2)-ztruthAll(1));
    coordsz1 = (ztruth+dzrand)/dz1+z01*ones(Nfits,1);
    coordsz2 = (ztruth+dzrand-zoffset)/dz2+z02*ones(Nfits,1);
    
    coords1 = [coordsxy1 coordsz1];% coordinates of the first channel
    coords2 = [coordsxy2 coordsz2];% coordinates of the second channel
    
    output1= simSplinePSF_call(Npixels,coeff1,Nphotons1,bg1,coords1);
    output2= simSplinePSF_call(Npixels,coeff2,Nphotons2,bg2,coords2);
    output1 = poissrnd(output1,Npixels,Npixels,Nfits);
    output2 = poissrnd(output2,Npixels,Npixels,Nfits);
    
    %fitting parameters
    iterations=50; % maximum fitting iterations
    fittype =2;% 2 for spline fit
    sCMOSvarmap = 0; %sCMOS noise map
    silent = 1;% display text, 1 if no text
    %data
    d_data(:,:,:,1) = output1;
    d_data(:,:,:,2) = output2;
    %spline coefficients
    coeff(:,:,:,:,1)=coeff1;
    coeff(:,:,:,:,2)=coeff2;
    
    %define shifts between channels
    noChannels = 2; % # of channels
    dT = zeros(5,noChannels,Nfits); 
    dxy=coordsxy1-coordsxy2;
    temp = reshape(dxy',[2 1,Nfits]);
    dT(1:2,2,:)=temp*-1;
    dT(3,2,:)=coordsz2-coordsz1; % parameter shifts between channels
    
    % define scaling between channels
    dS = repmat([1, 1 ;1, 1 ;1, 1;1, 1;1, 1],[1 1 Nfits]);
    
    % transformation between channels
    dTS = zeros(5,noChannels*2,Nfits);
    dTS(:,1:2,:)=dT;
    dTS(:,3:4,:)=dS;
    
    % linking parameters 
    shared_linkXYZ = repmat([1;1;1;0;0], [1 Nfits]); % only link XYZ
    shared_linkALL = repmat([1;1;1;1;1], [1 Nfits]); % link all parameters, XYZ,Photons and background
    
    
    
    % global fit link all parmeters
    [P_linkALL,CRLB_linkALL, LL1_linkALL] =  mleFit_LM_globalfit(d_data,fittype,shared_linkALL,iterations,coeff,dTS,sCMOSvarmap,silent,zstart);
    
    %global fit link xyz
    [P_linkXYZ,CRLB_linkXYZ, LL_linkXYZ] =  mleFit_LM_globalfit(d_data,fittype,shared_linkXYZ,iterations,coeff,dTS,sCMOSvarmap,silent,zstart);
    
    %individual fit
    [P_sep1,CRLB_sep1, LL_sep1] =  mleFit_LM(output1,5,50,coeff1,0,1,zstart(1,:));
    [P_sep2,CRLB_sep2, LL_sep2] =  mleFit_LM(output2,5,50,coeff2,0,1,zstart(1,:));
     %transform results back to channel 1 (reference coordiante system)
     xyback = zeros(Nfits,2);
    for i = 1:Nfits
        temp = tforminv*[P_sep2(i,1:2)';1];
        xyback(i,:) = temp(1:2);
    end
    Psep_Wcrkb=([P_sep1(:,1:2)./CRLB_sep1(:,1:2) P_sep1(:,5)./CRLB_sep1(:,5)]+[xyback./CRLB_sep2(:,1:2) P_sep2(:,5)./CRLB_sep2(:,5)])./(1./CRLB_sep1(:,[1 2 5])+1./CRLB_sep2(:,[1 2 5]));
    
    % ground truth parameters
    Theta(:,:,1) = [coordsxy1 coordsz1 Nphotons1*ones(Nfits,1) bg1*ones(Nfits,1)];
    Theta(:,:,2) = [coordsxy2 coordsz2 Nphotons2*ones(Nfits,1) bg2*ones(Nfits,1)];
   
    % CRLB with ground truth parameters
    [CRLBT_linkALL] =  calculate_CRLB_YL_multichannel_final(Nfits, coeff, Npixels, Theta,[1 1 1 1 1],dTS);
    [CRLBT_linkXYZ] =  calculate_CRLB_YL_multichannel_final(Nfits, coeff, Npixels, Theta,[1 1 1 0 0],dTS);
    
    
    
    %create filter for some bad localization data
    sigmazAll=mean(sqrt(CRLBT_linkALL(:,3)));
    temp = [coordsz1(1)-8*sigmazAll coordsz1(1)+8*sigmazAll];
    
    % RMSE_linkALL
    indP_linkALL = P_linkALL(:,3)>temp(1)&P_linkALL(:,3)<temp(2);
    fgood_linkALL(iz)=sum(indP_linkALL)/length(indP_linkALL);
    
    RMSEX_linkALLF(iz) = sqrt(mean((P_linkALL(indP_linkALL,1)-coordsxy1(indP_linkALL,1)).^2));
    RMSEY_linkALLF(iz) = sqrt(mean((P_linkALL(indP_linkALL,2)-coordsxy1(indP_linkALL,2)).^2));
    RMSEZ_linkALLF(iz) = sqrt(mean((P_linkALL(indP_linkALL,3)-coordsz1(indP_linkALL)).^2));
    RMSEP_linkALLF(iz) = sqrt(mean((P_linkALL(indP_linkALL,4)-Nphotons1).^2));
    
    msCRLBx_linkALLTF(iz) = mean(sqrt(CRLBT_linkALL(indP_linkALL,1)));
    msCRLBy_linkALLTF(iz) = mean(sqrt(CRLBT_linkALL(indP_linkALL,2)));
    msCRLBz_linkALLTF(iz) = mean(sqrt(CRLBT_linkALL(indP_linkALL,3)));
    msCRLBP_linkALLTF(iz) = mean(sqrt(CRLBT_linkALL(indP_linkALL,4)));
    
    %RMSE_linkXYZ
    indP_linkXYZ = P_linkXYZ(:,3)>temp(1)&P_linkXYZ(:,3)<temp(2);
    fgood_linkXYZ(iz)=sum(indP_linkXYZ)/length(indP_linkXYZ); 
    
    RMSEX_linkXYZF(iz) = sqrt(mean((P_linkXYZ(indP_linkXYZ,1)-coordsxy1(indP_linkXYZ,1)).^2));
    RMSEY_linkXYZF(iz) = sqrt(mean((P_linkXYZ(indP_linkXYZ,2)-coordsxy1(indP_linkXYZ,2)).^2));
    RMSEZ_linkXYZF(iz) = sqrt(mean((P_linkXYZ(indP_linkXYZ,3)-coordsz1(indP_linkXYZ)).^2));
    RMSEP1_linkXYZF(iz) = sqrt(mean((P_linkXYZ(indP_linkXYZ,4)-Nphotons1).^2));
    RMSEP2_linkXYZF(iz) = sqrt(mean((P_linkXYZ(indP_linkXYZ,5)-Nphotons2).^2));
    
    msCRLBx_linkXYZTF(iz) = mean(sqrt(CRLBT_linkXYZ(indP_linkXYZ,1)));
    msCRLBy_linkXYZTF(iz) = mean(sqrt(CRLBT_linkXYZ(indP_linkXYZ,2)));
    msCRLBz_linkXYZTF(iz) = mean(sqrt(CRLBT_linkXYZ(indP_linkXYZ,3)));
    msCRLBP1_linkXYZTF(iz) = mean(sqrt(CRLBT_linkXYZ(indP_linkXYZ,4)));
    msCRLBP2_linkXYZTF(iz) = mean(sqrt(CRLBT_linkXYZ(indP_linkXYZ,5)));
    
    %RMSE_individual Fit
    indP_sep1 = P_sep1(:,5)>temp(1)&P_sep1(:,5)<temp(2);
    indP_sep2 = P_sep2(:,5)>temp(1)&P_sep2(:,5)<temp(2);
    indP_sep = indP_sep1&indP_sep2;
   
    fgood_sep(iz)=sum(indP_sep)/length(indP_sep);
    
    RMSEX_WsepF(iz) = sqrt(mean((Psep_Wcrkb(indP_sep,1)-coordsxy1(indP_sep,1)).^2));
    RMSEY_WsepF(iz) = sqrt(mean((Psep_Wcrkb(indP_sep,2)-coordsxy1(indP_sep,2)).^2));
    RMSEZ_WsepF(iz) = sqrt(mean((Psep_Wcrkb(indP_sep,3)-coordsz1(indP_sep)).^2));
    RMSEP1_WsepF(iz) = sqrt(mean((P_sep1(indP_sep,3)-Nphotons1).^2));
    RMSEP2_WsepF(iz) = sqrt(mean((P_sep2(indP_sep,3)-Nphotons2).^2));
    
   

end


%% plot results

pixelSizeX = 100;
pixelSizeY = 100;
pixelSizeZ = 10;

figure(101)
hold off
plot(ztruthAll,msCRLBx_linkALLTF*pixelSizeX,'r')
hold on
plot(ztruthAll,RMSEX_linkALLF*pixelSizeX,'ro')
plot(ztruthAll,msCRLBy_linkALLTF*pixelSizeY,'g')
plot(ztruthAll,RMSEY_linkALLF*pixelSizeY,'go')
plot(ztruthAll,RMSEZ_linkALLF*pixelSizeZ,'bo')
plot(ztruthAll,msCRLBz_linkALLTF*pixelSizeZ,'b')
ylim([0 35])
xlim([-500 500])
xlabel('z (nm)')
ylabel('localization error (nm)')
title('global fit all linked')


figure(102)
hold off
plot(ztruthAll,msCRLBx_linkXYZTF*pixelSizeX,'r')
hold on
plot(ztruthAll,RMSEX_linkXYZF*pixelSizeX,'ro')
plot(ztruthAll,msCRLBy_linkXYZTF*pixelSizeY,'g')
plot(ztruthAll,RMSEY_linkXYZF*pixelSizeY,'go')
plot(ztruthAll,RMSEZ_linkXYZF*pixelSizeZ,'bo')
plot(ztruthAll,msCRLBz_linkXYZTF*pixelSizeZ,'b')
ylim([0 35])
xlim([-500 500])
xlabel('z (nm)')
ylabel('localization error (nm)')
title('global fit xyz linked')



figure(103)
hold off
plot(ztruthAll,msCRLBx_linkALLTF*pixelSizeX,'r')
hold on
plot(ztruthAll,RMSEX_WsepF*pixelSizeX,'ro')
plot(ztruthAll,msCRLBy_linkALLTF*pixelSizeY,'g')
plot(ztruthAll,RMSEY_WsepF*pixelSizeY,'go')
plot(ztruthAll,RMSEZ_WsepF*pixelSizeZ,'bo')
plot(ztruthAll,msCRLBz_linkALLTF*pixelSizeZ,'b')
ylim([0 35])
xlim([-500 500])
xlabel('z (nm)')
ylabel('localization error (nm)')
title('individual fits')


figure(104)
plot(ztruthAll, fgood_linkALL,ztruthAll,fgood_linkXYZ,ztruthAll,fgood_sep);
xlabel('z (nm)')
ylabel('fraction of good fits')
legend('all linked','xyz linked','individual fit')






















