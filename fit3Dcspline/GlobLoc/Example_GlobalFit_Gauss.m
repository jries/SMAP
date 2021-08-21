%% test Gauss multi channel fit


Nfits = 1000;
Nphotons1 =2500;
Nphotons2 =2500;
Npixels = 15;
bg1 = 20;
bg2 = 20;
zoffset = 0;%nm

% ztruthAll = -600:50:600;
% numz=length(ztruthAll);

% zstartnm = 200;
% zstart=repmat([-2 -1 1 2]*zstartnm/dz1,[Nfits 1]);
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

PSFsigma1=1.2;
PSFsigma2=2;%PSF sigma in pixels
fittype=1;


coordsxy1 = Npixels/2 -1 +2*rand([Nfits 2]);
coordsxy2 = zeros(Nfits,2);

for i = 1:Nfits
    temp = tformF*[coordsxy1(i,:)';1];
    coordsxy2(i,:) = temp(1:2);
end


%   Generate a stack of images
% coords=Npixels/2-1+rand([Nfits 2])+[0*ones(Nfits,1) zeros(Nfits,1)];
[output1] = finitegausspsf(Npixels,PSFsigma1,Nphotons1,bg1,coordsxy1);
[output2] = finitegausspsf(Npixels,PSFsigma2,Nphotons2,bg2,coordsxy2);
output1 = poissrnd(output1,Npixels,Npixels,Nfits);
output2 = poissrnd(output2,Npixels,Npixels,Nfits);
d_data(:,:,:,1) = output1;
d_data(:,:,:,2) = output2;

noChannels = 2; % # of channels
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
iterations = 50;
sCMOSvarmap = 0;
silent = 0;

shared_linkALL = repmat([1;1;1;1;1], [1 Nfits]); % link all parameters, XYZ,Photons and background




 [P_linkALL,CRLB_linkALL, LL1_linkALL] =  mleFit_LM_globalfit(d_data,fittype,shared_linkALL,iterations,single(1),dTS,sCMOSvarmap,silent);







%% Gauss

Nfits=500    %number of images to fit
bg=10;           %background fluorescence in photons/pixel/frame
Nphotons=2500;   %expected photons/frame
% Npixels=9;      %linear size of fit region in pixels. 
PSFsigma1=1.2; 
PSFsigma2=1.5;%PSF sigma in pixels
fittype=2;

n=1;
for Npixels = 7:19
    %   Generate a stack of images
    coords=Npixels/2-1+rand([Nfits 2])+[0*ones(Nfits,1) zeros(Nfits,1)];
    [out] = finitegausspsf(Npixels,PSFsigma1,Nphotons,bg,coords);
    
    %   Corrupt with Poisson noise
    data = poissrnd(out,Npixels,Npixels,Nfits);
    %data = poissrnd(out); %requires statistics toolbox
    %   Can look at data (DipImage)
    %dipshow(permute(data,[2 1 3])); %dipimage permuted 1st two dimensions
    
    coords1 = coords;
    data1 = data;
    
    % coords=Npixels/2-1+rand([Nfits 2])+[0*ones(Nfits,1) zeros(Nfits,1)];
    coords = coords1+1;
    [out] = finitegausspsf(Npixels,PSFsigma2,Nphotons,bg,coords);
    
    %   Corrupt with Poisson noise
%     data=single(noise(out,'poisson',1)); %requires DipImage
    data = poissrnd(out,Npixels,Npixels,Nfits);
% output2 = poissrnd(output2,Npixels,Npixels,Nfits);
    coords2 = coords;
    data2 = data;
    
    
    d_data(:,:,:,1) = data1;
    d_data(:,:,:,2) = data2;
    
    
    noChannels = 2;
    iterations=50;
    dT = zeros(5,noChannels,Nfits);
    dxy=coords1-coords2;
    temp = reshape(dxy',[2 1,Nfits]);
    dT(1:2,2,:)=temp*(-1);
    % dT(5,2,:)=0.8;
    
    shared = [1;1;1;1;0];
    sharedA = uint32(repmat(shared,[1 Nfits]));
    fittype = uint32(2);
    initSigma = single(1.2);
    dT = single(dT);
    d_data(:,:,:,1) = data1;
    d_data(:,:,:,2) = data2;
    

    tic
%     [PG,CRLBG, LLG] =  GPUmleFit_LM_MultiChannel(d_data,fittype,(sharedA),50,initSigma,(dT));
    [PGC,CRLBGC, LLGC] =  CPUmleFit_LM_MultiChannel(single(d_data),fittype,(sharedA),50,initSigma,(dT));

    tGauss=toc;
    
    GPUmleFit_LM_MultiChannel_Gauss(n,1)=Nfits/tGauss;
    GPUmleFit_LM_MultiChannel_Gauss(n,2)=Npixels;
    
%     CPUmleFit_LM_MultiChannel_Gauss(n,1)=Nfits/tGauss;
%     CPUmleFit_LM_MultiChannel_Gauss(n,2)=Npixels;
    
    
    clear d_data
    n = n+1
    
end











