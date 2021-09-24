% Copyright (c) 2021 Li Lab, Southern University of Science and Technology, Shenzhen & Ries Lab, European Molecular Biology Laboratory, Heidelberg.
% author: Yiming Li
% email: liym2019@sustech.edu.cn
% date: 2021.08.27
% Tested with CUDA 11.3 (Express installation) and Matlab 2019a
%% Gauss
clearvars
Nfits=10000    %number of images to fit
bg=20;           %background fluorescence in photons/pixel/frame
Nphotons1=2500;   %expected photons/frame
Nphotons2 =2500;
bg1 = 20;
bg2 = 20;
% Npixels=9;      %linear size of fit region in pixels. 
PSFsigma1=1.2; 
PSFsigma2=1.5;%PSF sigma in pixels
fittype=1;%1 for Gauss fit, 2 for spline fit

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

noChannels = 2; %number of channels
iterations=50;  %iteration time
sCMOSvarmap = 0;
silent = 0;
initSigma = 1.5;
n=1;
for Npixels = 7:19
    %   Generate a stack of images
    coordsxy1 = Npixels/2 -1 +2*rand([Nfits 2]);
    coordsxy2 = zeros(Nfits,2);
    for i = 1:Nfits
        temp = tformF*[coordsxy1(i,:)';1];
        coordsxy2(i,:) = temp(1:2);
    end
    
    
    [output1] = finitegausspsf(Npixels,PSFsigma1,Nphotons1,bg1,coordsxy1);
    [output2] = finitegausspsf(Npixels,PSFsigma2,Nphotons2,bg2,coordsxy2);
    output1 = poissrnd(output1,Npixels,Npixels,Nfits);
    output2 = poissrnd(output2,Npixels,Npixels,Nfits);
    d_data(:,:,:,1) = output1;
    d_data(:,:,:,2) = output2;
    
    
    
    dT = zeros(5,noChannels,Nfits);
    dxy=coordsxy1-coordsxy2;
    temp = reshape(dxy',[2 1,Nfits]);
    dT(1:2,2,:)=temp*-1;
    dT(5,2,:)=PSFsigma2-PSFsigma1; % parameter shifts between channels
    
    % define scaling between channels
    dS = repmat([1, 1 ;1, 1 ;1, 1;1, 1;1, 1],[1 1 Nfits]);
    
    % transformation between channels
    dTS = zeros(5,noChannels*2,Nfits);
    dTS(:,1:2,:)=dT;
    dTS(:,3:4,:)=dS;

    
    shared = repmat([1;1;1;1;1], [1 Nfits]); % link all parameters, XYZ,Photons and background

    
   
    d_data = single(d_data);
    fittype = int32(fittype);
    shared = int32(shared);
    iterations = int32(iterations);
    initSigma = single(initSigma);
    dTS = single(dTS);

    tic
%      [P,CRLB, LL1] =  mleFit_LM_globalfit(d_data,fittype,shared,iterations,initSigma,dTS,sCMOSvarmap,silent);
     
     [P,CRLB, LL] =  GPUmleFit_LM_MultiChannel(d_data,fittype,shared,iterations,initSigma,dTS,0,1);
%      [P,CRLB, LL] =  CPUmleFit_LM_MultiChannel(d_data,fittype,shared,iterations,initSigma,dTS,0,1);
     
    tGauss=toc;
    
    GPUmleFit_LM_MultiChannel_Gauss(n,1)=Nfits/tGauss;
    GPUmleFit_LM_MultiChannel_Gauss(n,2)=Npixels;
    disp(['Actual Gaussian fits per second in GPU:' num2str(Nfits/tGauss) ' for ROI ' num2str(Npixels) ' pixels'])
%     CPUmleFit_LM_MultiChannel_Gauss(n,1)=Nfits/tGauss;
%     CPUmleFit_LM_MultiChannel_Gauss(n,2)=Npixels;
%     disp(['Actual Gaussian fits per second in CPU:' num2str(Nfits/tGauss) ' for ROI ' num2str(Npixels) ' pixels'])
    
    
    clear d_data
    n = n+1;
    
end











