% Copyright (c) 2021 Li Lab, Southern University of Science and Technology, Shenzhen & Ries Lab, European Molecular Biology Laboratory, Heidelberg.
% author: Yiming Li
% email: liym2019@sustech.edu.cn
% date: 2021.08.27
% Download the calbration data for the followling link first: https://www.embl.de/download/ries/globLoc/
%000_AstigBeads_LP665_m1_1um_20nm_50ms_conv_640_50Percent_Empty_685-70_676-37_singleMode_Z-stack_1_MMStack_Pos0.ome_3dcal.mat
% Tested with CUDA 11.3 (Express installation) and Matlab 2019a
%% Load and make cspline PSF
clearvars
cal=load('000_AstigBeads_LP665_m1_1um_20nm_50ms_conv_640_50Percent_Empty_685-70_676-37_singleMode_Z-stack_1_MMStack_Pos0.ome_3dcal.mat');
PSF1 = cal.SXY(1).PSF{1};
rangex = 4:30;
rangey = 4:30;
rangez = 2:100;
z = 51;
corrPSFr = PSF1;
centpsfr=PSF1(rangex,rangey,z-1:z+1); %cut out rim from shift
minPSFr=min(centpsfr(:),[],'omitnan');
corrPSFnr=corrPSFr-minPSFr;
intglobalr=mean(sum(sum(corrPSFnr(rangex,rangey,z-1:z+1),1,'omitnan'),2,'omitnan'),'omitnan');
corrPSFnr=corrPSFnr/intglobalr;
corrPSFnr(isnan(corrPSFnr))=0;
corrPSFnr(corrPSFnr<0)=0;
corrPSFsr=corrPSFnr(rangex,rangey,rangez);
bsplinePSF1=bsarray(double(corrPSFsr),'lambda',[0 0 1]);

% zhd=1:1:b3_0r.dataSize(3);
% dxxhd=1;
[XX,YY,ZZ]=meshgrid(1:bsplinePSF1.dataSize(1),1:bsplinePSF1.dataSize(2),1:bsplinePSF1.dataSize(3));
 PSFsmooth1 = interp3_0(bsplinePSF1,XX,YY,ZZ,0);

coeff1 = single(Spline3D_interp(PSFsmooth1));
% coeff1 = cal.SXY(1).cspline.coeff{1};
dz1=cal.SXY_g.cspline.dz;
z01=cal.SXY_g.cspline.z0-8;%shift the stack to focus


PSF2 = cal.SXY(2).PSF{1};
rangex = 4:30;
rangey = 4:30;
rangez = 2:100;
z = 51;
corrPSFr = PSF2;
centpsfr=PSF2(rangex,rangey,z-1:z+1); %cut out rim from shift
minPSFr=min(centpsfr(:),[],'omitnan');
corrPSFnr=corrPSFr-minPSFr;
intglobalr=mean(sum(sum(corrPSFnr(rangex,rangey,z-1:z+1),1,'omitnan'),2,'omitnan'),'omitnan');
corrPSFnr=corrPSFnr/intglobalr;
corrPSFnr(isnan(corrPSFnr))=0;
corrPSFnr(corrPSFnr<0)=0;
corrPSFsr=corrPSFnr(rangex,rangey,rangez);
bsplinePSF2=bsarray(double(corrPSFsr),'lambda',[0 0 1]);


[XX,YY,ZZ]=meshgrid(1:bsplinePSF2.dataSize(1),1:bsplinePSF2.dataSize(2),1:bsplinePSF2.dataSize(3));
 PSFsmooth2 = interp3_0(bsplinePSF2,XX,YY,ZZ,0);

coeff2 = single(Spline3D_interp(PSFsmooth2));
dz2=cal.SXY_g.cspline.dz;
z02=cal.SXY_g.cspline.z0-8; %shift the stack to focus

%% Parameters
%load photon distributions of 4 dyes
distribution = 'F_Lxyz_Nup96SNAPAF647_ELYSCF660C_NUP62DY634_WGA_CF680_100p_10ms_EM100_2_Pos0_driftc_photons';
load([distribution '.mat'])
Nfits = 5000;
Nphotons1up = rand_arb_cjw(Nfits, out(:,[1 2])); 
Nphotons1down = 1466/3798*Nphotons1up;%DY634
Nphotons2up = rand_arb_cjw(Nfits, out(:,[1 3])); 
Nphotons2down = 1987/9350*Nphotons2up;%AF647
Nphotons3up = rand_arb_cjw(Nfits, out(:,[1 4])); 
Nphotons3down = 438/6423*Nphotons3up;%CF660C
Nphotons4up = rand_arb_cjw(Nfits, out(:,[1 5])); 
Nphotons4down = 126/7061*Nphotons4up;%CF680

%define affine transformation matrix
theta = 10*pi/180;%rotation
sx = 1.01;%scalex
sy = 1.02;%scaley
tx = -.5;%shfitx
ty = -0.2;%shfity
tformR = [cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0;0 0 1];
tformS = [sx 0 0;0 sy 0;0 0 1];
tformT = [1 0 tx; 0 1 ty; 0 0 1];
tform = tformR*tformS*tformT;
tforminv=inv(tform);


% simulate data
Npixels = 15;
bg1 = 20;
bg2 = 20;
bg3 = 20;
bg4 = 20;
zoffset = 0;% focus difference between two channels nm
coordsxy1 = Npixels/2 -1 +2*rand([Nfits 2]);
coordsxy2 = zeros(Nfits,2);
for i = 1:Nfits
    temp = tform*[coordsxy1(i,:)';1];
    coordsxy2(i,:) = temp(1:2);
end
ztruth = -600+1200*rand(Nfits,1);
coordsz1 = ztruth/dz1+z01*ones(Nfits,1);
coordsz2 = (ztruth-zoffset)/dz2+z02*ones(Nfits,1);

coords1 = [coordsxy1 coordsz1];
coords2 = [coordsxy2 coordsz2];

output1up = simSplinePSF_call(Npixels,coeff1,Nphotons1up,bg1,coords1);
output1down = simSplinePSF_call(Npixels,coeff2,Nphotons1down,bg1,coords2);

output2up = simSplinePSF_call(Npixels,coeff1,Nphotons2up,bg2,coords1);
output2down = simSplinePSF_call(Npixels,coeff2,Nphotons2down,bg2,coords2);

output3up = simSplinePSF_call(Npixels,coeff1,Nphotons3up,bg3,coords1);
output3down = simSplinePSF_call(Npixels,coeff2,Nphotons3down,bg3,coords2);

output4up = simSplinePSF_call(Npixels,coeff1,Nphotons4up,bg4,coords1);
output4down = simSplinePSF_call(Npixels,coeff2,Nphotons4down,bg4,coords2);

output1up = poissrnd(output1up,Npixels,Npixels,Nfits);
output1down = poissrnd(output1down,Npixels,Npixels,Nfits);

output2up = poissrnd(output2up,Npixels,Npixels,Nfits);
output2down = poissrnd(output2down,Npixels,Npixels,Nfits);

output3up = poissrnd(output3up,Npixels,Npixels,Nfits);
output3down = poissrnd(output3down,Npixels,Npixels,Nfits);

output4up = poissrnd(output4up,Npixels,Npixels,Nfits);
output4down = poissrnd(output4down,Npixels,Npixels,Nfits);

d_data1(:,:,:,1) = output1up;
d_data1(:,:,:,2) = output1down;

d_data2(:,:,:,1) = output2up;
d_data2(:,:,:,2) = output2down;

d_data3(:,:,:,1) = output3up;
d_data3(:,:,:,2) = output3down;

d_data4(:,:,:,1) = output4up;
d_data4(:,:,:,2) = output4down;

coeff(:,:,:,:,1)=coeff1;
coeff(:,:,:,:,2)=coeff2;


%fitting parameters
iterations=50; % maximum fitting iterations
fittype =2;% 2 for spline fit
sCMOSvarmap = 0; %sCMOS noise map
silent = 0;% display text, 1 if no text

%define shifts between channels
noChannels = 2; % # of channels
dT = zeros(5,noChannels,Nfits);
dxy=coordsxy1-coordsxy2;
temp = reshape(dxy',[2 1,Nfits]);
dT(1:2,2,:)=temp*-1;
dT(3,2,:)=coordsz2-coordsz1; % parameter shifts between channels

% define scaling between channels
dS1 = repmat([1, 1 ;1, 1 ;1, 1;1, 1466/3798;1, 1],[1 1 Nfits]);
dS2 = repmat([1, 1 ;1, 1 ;1, 1;1, 1987/9350;1, 1],[1 1 Nfits]);
dS3 = repmat([1, 1 ;1, 1 ;1, 1;1, 438/6423;1, 1],[1 1 Nfits]);
dS4 = repmat([1, 1 ;1, 1 ;1, 1;1, 126/7061;1, 1],[1 1 Nfits]);

% transformation between channels
dTS1 = zeros(5,noChannels*2,Nfits);
dTS2 = zeros(5,noChannels*2,Nfits);
dTS3 = zeros(5,noChannels*2,Nfits);
dTS4 = zeros(5,noChannels*2,Nfits);
dTS1(:,1:2,:)=dT;
dTS1(:,3:4,:)=dS1;
dTS2(:,1:2,:)=dT;
dTS2(:,3:4,:)=dS2;
dTS3(:,1:2,:)=dT;
dTS3(:,3:4,:)=dS3;
dTS4(:,1:2,:)=dT;
dTS4(:,3:4,:)=dS4;

% linking parameters
shared_linkXYZ = repmat([1;1;1;0;0], [1 Nfits]); % only link XYZ
shared_linkXYZP = repmat([1;1;1;1;0], [1 Nfits]); % link all parameters, XYZ,Photons and background

zstartnm = 200;
zstart=repmat([-2 -1 1 2]*zstartnm/dz1,[Nfits 1]);

%Ratiometric Fit**************************************************************************************************************************************************************
%fit dye1 with different ratios
[P1_LinkP1,CRLB1_LinkP1, LL1_LinkP1] =  mleFit_LM_globalfit(d_data1,fittype,shared_linkXYZP,iterations,coeff,dTS1,sCMOSvarmap,silent,zstart);
[P1_LinkP2,CRLB1_LinkP2, LL1_LinkP2] =  mleFit_LM_globalfit(d_data1,fittype,shared_linkXYZP,iterations,coeff,dTS2,sCMOSvarmap,silent,zstart);
[P1_LinkP3,CRLB1_LinkP3, LL1_LinkP3] =  mleFit_LM_globalfit(d_data1,fittype,shared_linkXYZP,iterations,coeff,dTS3,sCMOSvarmap,silent,zstart);
[P1_LinkP4,CRLB1_LinkP4, LL1_LinkP4] =  mleFit_LM_globalfit(d_data1,fittype,shared_linkXYZP,iterations,coeff,dTS4,sCMOSvarmap,silent,zstart);

PALL1(:,:,1)= P1_LinkP1;
PALL1(:,:,2)= P1_LinkP2;
PALL1(:,:,3)= P1_LinkP3;
PALL1(:,:,4)= P1_LinkP4;
P_linkP1 = zeros(size(P1_LinkP1));

CRLBALL1(:,:,1)= CRLB1_LinkP1;
CRLBALL1(:,:,2)= CRLB1_LinkP2;
CRLBALL1(:,:,3)= CRLB1_LinkP3;
CRLBALL1(:,:,4)= CRLB1_LinkP4;
CRLB_linkP1 = zeros(size(CRLB1_LinkP1));

LLALL1 = [LL1_LinkP1 LL1_LinkP2 LL1_LinkP3 LL1_LinkP4];

LL_linkP1 = zeros(size(LL1_LinkP1));

[M1,I1]=max(LLALL1,[],2);
for i = 1:Nfits
    P_linkP1(i,:) = squeeze(PALL1(i,:,I1(i)));
    CRLB_linkP1(i,:) = squeeze(CRLBALL1(i,:,I1(i)));
    LL_linkP1(i)=LLALL1(i,I1(i));
end


%fit dye2 with different ratios
[P2_LinkP1,CRLB2_LinkP1, LL2_LinkP1] =  mleFit_LM_globalfit(d_data2,fittype,shared_linkXYZP,iterations,coeff,dTS1,sCMOSvarmap,silent,zstart);
[P2_LinkP2,CRLB2_LinkP2, LL2_LinkP2] =  mleFit_LM_globalfit(d_data2,fittype,shared_linkXYZP,iterations,coeff,dTS2,sCMOSvarmap,silent,zstart);
[P2_LinkP3,CRLB2_LinkP3, LL2_LinkP3] =  mleFit_LM_globalfit(d_data2,fittype,shared_linkXYZP,iterations,coeff,dTS3,sCMOSvarmap,silent,zstart);
[P2_LinkP4,CRLB2_LinkP4, LL2_LinkP4] =  mleFit_LM_globalfit(d_data2,fittype,shared_linkXYZP,iterations,coeff,dTS4,sCMOSvarmap,silent,zstart);

PALL2(:,:,1)= P2_LinkP1;
PALL2(:,:,2)= P2_LinkP2;
PALL2(:,:,3)= P2_LinkP3;
PALL2(:,:,4)= P2_LinkP4;
P_linkP2 = zeros(size(P2_LinkP1));

CRLBALL2(:,:,1)= CRLB2_LinkP1;
CRLBALL2(:,:,2)= CRLB2_LinkP2;
CRLBALL2(:,:,3)= CRLB2_LinkP3;
CRLBALL2(:,:,4)= CRLB2_LinkP4;
CRLB_linkP2 = zeros(size(CRLB2_LinkP1));

LLALL2 = [LL2_LinkP1 LL2_LinkP2 LL2_LinkP3 LL2_LinkP4];

LL_linkP2 = zeros(size(LL2_LinkP1));

[M2,I2]=max(LLALL2,[],2);
for i = 1:Nfits
    P_linkP2(i,:) = squeeze(PALL2(i,:,I2(i)));
    CRLB_linkP2(i,:) = squeeze(CRLBALL2(i,:,I2(i)));
    LL_linkP2(i)=LLALL2(i,I2(i));
    
end

%fit dye3 with different ratios
[P3_LinkP1,CRLB3_LinkP1, LL3_LinkP1] =  mleFit_LM_globalfit(d_data3,fittype,shared_linkXYZP,iterations,coeff,dTS1,sCMOSvarmap,silent,zstart);
[P3_LinkP2,CRLB3_LinkP2, LL3_LinkP2] =  mleFit_LM_globalfit(d_data3,fittype,shared_linkXYZP,iterations,coeff,dTS2,sCMOSvarmap,silent,zstart);
[P3_LinkP3,CRLB3_LinkP3, LL3_LinkP3] =  mleFit_LM_globalfit(d_data3,fittype,shared_linkXYZP,iterations,coeff,dTS3,sCMOSvarmap,silent,zstart);
[P3_LinkP4,CRLB3_LinkP4, LL3_LinkP4] =  mleFit_LM_globalfit(d_data3,fittype,shared_linkXYZP,iterations,coeff,dTS4,sCMOSvarmap,silent,zstart);

PALL3(:,:,1)= P3_LinkP1;
PALL3(:,:,2)= P3_LinkP2;
PALL3(:,:,3)= P3_LinkP3;
PALL3(:,:,4)= P3_LinkP4;
P_linkP3 = zeros(size(P3_LinkP1));

CRLBALL3(:,:,1)= CRLB3_LinkP1;
CRLBALL3(:,:,2)= CRLB3_LinkP2;
CRLBALL3(:,:,3)= CRLB3_LinkP3;
CRLBALL3(:,:,4)= CRLB3_LinkP4;
CRLB_linkP3 = zeros(size(CRLB3_LinkP1));

LLALL3 = [LL3_LinkP1 LL3_LinkP2 LL3_LinkP3 LL3_LinkP4];

LL_linkP3 = zeros(size(LL3_LinkP1));

[M3,I3]=max(LLALL3,[],2);
for i = 1:Nfits
    P_linkP3(i,:) = squeeze(PALL3(i,:,I3(i)));
    CRLB_linkP3(i,:) = squeeze(CRLBALL3(i,:,I3(i)));
    LL_linkP3(i)=LLALL3(i,I3(i));
    
end

%fit dye4 with different ratios
[P4_LinkP1,CRLB4_LinkP1, LL4_LinkP1] =  mleFit_LM_globalfit(d_data4,fittype,shared_linkXYZP,iterations,coeff,dTS1,sCMOSvarmap,silent,zstart);
[P4_LinkP2,CRLB4_LinkP2, LL4_LinkP2] =  mleFit_LM_globalfit(d_data4,fittype,shared_linkXYZP,iterations,coeff,dTS2,sCMOSvarmap,silent,zstart);
[P4_LinkP3,CRLB4_LinkP3, LL4_LinkP3] =  mleFit_LM_globalfit(d_data4,fittype,shared_linkXYZP,iterations,coeff,dTS3,sCMOSvarmap,silent,zstart);
[P4_LinkP4,CRLB4_LinkP4, LL4_LinkP4] =  mleFit_LM_globalfit(d_data4,fittype,shared_linkXYZP,iterations,coeff,dTS4,sCMOSvarmap,silent,zstart);

PALL4(:,:,1)= P4_LinkP1;
PALL4(:,:,2)= P4_LinkP2;
PALL4(:,:,3)= P4_LinkP3;
PALL4(:,:,4)= P4_LinkP4;
P_linkP4 = zeros(size(P4_LinkP1));

CRLBALL4(:,:,1)= CRLB4_LinkP1;
CRLBALL4(:,:,2)= CRLB4_LinkP2;
CRLBALL4(:,:,3)= CRLB4_LinkP3;
CRLBALL4(:,:,4)= CRLB4_LinkP4;
CRLB_linkP4 = zeros(size(CRLB4_LinkP1));

LLALL4 = [LL4_LinkP1 LL4_LinkP2 LL4_LinkP3 LL4_LinkP4];

LL_linkP4 = zeros(size(LL4_LinkP1));

[M4,I4]=max(LLALL4,[],2);
for i = 1:Nfits
    P_linkP4(i,:) = squeeze(PALL4(i,:,I4(i)));
    CRLB_linkP4(i,:) = squeeze(CRLBALL4(i,:,I4(i)));
    LL_linkP4(i)=LLALL4(i,I4(i));    
end

% separation of ratiometric fit
ratioThreshold =0.999;% if the ratio between the best and second best LL are within ratioThreshold, the fit is discarded

%calculate crosstalk of each dyes
maxLLALL1=max(LLALL1,[],2);
maxLLALL1 = repmat(maxLLALL1,1,4);
ratioLLALL1 = maxLLALL1./LLALL1;
ratioLLALL11 = sort(ratioLLALL1,2);
ind1 = ratioLLALL11(:,3)<ratioThreshold;
I1F = I1(ind1);
cross11_linkPF = size(I1F(I1F==1),1)/size(I1F,1);
cross12_linkPF = size(I1F(I1F==2),1)/size(I1F,1);
cross13_linkPF = size(I1F(I1F==3),1)/size(I1F,1);
cross14_linkPF = size(I1F(I1F==4),1)/size(I1F,1);

maxLLALL2=max(LLALL2,[],2);
maxLLALL2 = repmat(maxLLALL2,1,4);
ratioLLALL2 = maxLLALL2./LLALL2;
ratioLLALL21 = sort(ratioLLALL2,2);
ind2 = ratioLLALL21(:,3)<ratioThreshold;
I2F = I2(ind2);
cross21_linkPF = size(I2F(I2F==1),1)/size(I2F,1);
cross22_linkPF = size(I2F(I2F==2),1)/size(I2F,1);
cross23_linkPF = size(I2F(I2F==3),1)/size(I2F,1);
cross24_linkPF = size(I2F(I2F==4),1)/size(I2F,1);

maxLLALL3=max(LLALL3,[],2);
maxLLALL3 = repmat(maxLLALL3,1,4);
ratioLLALL3 = maxLLALL3./LLALL3;
ratioLLALL31 = sort(ratioLLALL3,2);
ind3 = ratioLLALL31(:,3)<ratioThreshold;
I3F = I3(ind3);
cross31_linkPF = size(I3F(I3F==1),1)/size(I3F,1);
cross32_linkPF = size(I3F(I3F==2),1)/size(I3F,1);
cross33_linkPF = size(I3F(I3F==3),1)/size(I3F,1);
cross34_linkPF = size(I3F(I3F==4),1)/size(I3F,1);

maxLLALL4=max(LLALL4,[],2);
maxLLALL4 = repmat(maxLLALL4,1,4);
ratioLLALL4 = maxLLALL4./LLALL4;
ratioLLALL41 = sort(ratioLLALL4,2);
ind4 = ratioLLALL41(:,3)<ratioThreshold;
I4F = I4(ind4);
cross41_linkPF = size(I4F(I4F==1),1)/size(I4F,1);
cross42_linkPF = size(I4F(I4F==2),1)/size(I4F,1);
cross43_linkPF = size(I4F(I4F==3),1)/size(I4F,1);
cross44_linkPF = size(I4F(I4F==4),1)/size(I4F,1);

%percentage of fit kept
goodnessI1_linkP = sum(ind1)/size(ind1,1);
goodnessI2_linkP = sum(ind2)/size(ind2,1);
goodnessI3_linkP = sum(ind3)/size(ind3,1);
goodnessI4_linkP = sum(ind4)/size(ind4,1);
goodnessALL_linkXYZP =[goodnessI1_linkP;goodnessI2_linkP;goodnessI3_linkP;goodnessI4_linkP];


% fit with only linkXYZ**************************************************************************************************************************************************************
% dye1
[P1,CRLB1, LL1] =  mleFit_LM_globalfit(d_data1,fittype,shared_linkXYZ,iterations,coeff,dTS1,sCMOSvarmap,silent,zstart);
% dye2
[P2,CRLB2, LL2] =  mleFit_LM_globalfit(d_data2,fittype,shared_linkXYZ,iterations,coeff,dTS2,sCMOSvarmap,silent,zstart);
% dye3 
[P3,CRLB3, LL3] =  mleFit_LM_globalfit(d_data3,fittype,shared_linkXYZ,iterations,coeff,dTS3,sCMOSvarmap,silent,zstart);
% dye4 
[P4,CRLB4, LL4] =  mleFit_LM_globalfit(d_data4,fittype,shared_linkXYZ,iterations,coeff,dTS4,sCMOSvarmap,silent,zstart);

% calculate ratio
ratio1_linkXYZ = (P1(:,4)-P1(:,5))./(P1(:,4)+P1(:,5));
ratio2_linkXYZ = (P2(:,4)-P2(:,5))./(P2(:,4)+P2(:,5));
ratio3_linkXYZ = (P3(:,4)-P3(:,5))./(P3(:,4)+P3(:,5));
ratio4_linkXYZ = (P4(:,4)-P4(:,5))./(P4(:,4)+P4(:,5));
ratio_LinkXYZ_ALL = [ratio1_linkXYZ;ratio2_linkXYZ;ratio3_linkXYZ;ratio4_linkXYZ];

%define threshold
threshold1 = 0.5457;
threshold2 = 0.7668;
threshold3 = 0.9116;

IntervalT = 0.01;
threshold1 = [0.5457-IntervalT 0.5457+IntervalT];
threshold2 = [0.7668-IntervalT 0.7668+IntervalT];
threshold3 = [0.9116-IntervalT 0.9116+IntervalT];

% separaration of linkXYZ
color1_linkXYZ = ratio1_linkXYZ(ratio1_linkXYZ<=threshold1(1)|(ratio1_linkXYZ>threshold1(2)&ratio1_linkXYZ<=threshold2(1))|(ratio1_linkXYZ>threshold2(2)&ratio1_linkXYZ<=threshold3(1)|ratio1_linkXYZ>threshold3(2)));
cross11_linkXYZ = size(ratio1_linkXYZ(ratio1_linkXYZ<=threshold1(1)),1)/size(color1_linkXYZ,1);
cross12_linkXYZ=size(ratio1_linkXYZ(ratio1_linkXYZ>threshold1(2)&ratio1_linkXYZ<=threshold2(1)),1)/size(color1_linkXYZ,1);
cross13_linkXYZ=size(ratio1_linkXYZ(ratio1_linkXYZ>threshold2(2)&ratio1_linkXYZ<=threshold3(1)),1)/size(color1_linkXYZ,1);
cross14_linkXYZ=size(ratio1_linkXYZ(ratio1_linkXYZ>threshold3(2)),1)/size(color1_linkXYZ,1);

color2_linkXYZ = ratio2_linkXYZ(ratio2_linkXYZ<=threshold1(1)|(ratio2_linkXYZ>threshold1(2)&ratio2_linkXYZ<=threshold2(1))|(ratio2_linkXYZ>threshold2(2)&ratio2_linkXYZ<=threshold3(1)|ratio2_linkXYZ>threshold3(2)));
cross21_linkXYZ = size(ratio2_linkXYZ(ratio2_linkXYZ<=threshold1(1)),1)/size(color2_linkXYZ,1);
cross22_linkXYZ=size(ratio2_linkXYZ(ratio2_linkXYZ>threshold1(2)&ratio2_linkXYZ<=threshold2(1)),1)/size(color2_linkXYZ,1);
cross23_linkXYZ=size(ratio2_linkXYZ(ratio2_linkXYZ>threshold2(2)&ratio2_linkXYZ<=threshold3(1)),1)/size(color2_linkXYZ,1);
cross24_linkXYZ=size(ratio2_linkXYZ(ratio2_linkXYZ>threshold3(2)),1)/size(color2_linkXYZ,1);

color3_linkXYZ = ratio3_linkXYZ(ratio3_linkXYZ<=threshold1(1)|(ratio3_linkXYZ>threshold1(2)&ratio3_linkXYZ<=threshold2(1))|(ratio3_linkXYZ>threshold2(2)&ratio3_linkXYZ<=threshold3(1)|ratio3_linkXYZ>threshold3(2)));
cross31_linkXYZ = size(ratio3_linkXYZ(ratio3_linkXYZ<=threshold1(1)),1)/size(color3_linkXYZ,1);
cross32_linkXYZ=size(ratio3_linkXYZ(ratio3_linkXYZ>threshold1(2)&ratio3_linkXYZ<=threshold2(1)),1)/size(color3_linkXYZ,1);
cross33_linkXYZ=size(ratio3_linkXYZ(ratio3_linkXYZ>threshold2(2)&ratio3_linkXYZ<=threshold3(1)),1)/size(color3_linkXYZ,1);
cross34_linkXYZ=size(ratio3_linkXYZ(ratio3_linkXYZ>threshold3(2)),1)/size(color3_linkXYZ,1);

color4_linkXYZ = ratio4_linkXYZ(ratio4_linkXYZ<=threshold1(1)|(ratio4_linkXYZ>threshold1(2)&ratio4_linkXYZ<=threshold2(1))|(ratio4_linkXYZ>threshold2(2)&ratio4_linkXYZ<=threshold3(1)|ratio4_linkXYZ>threshold3(2)));
cross41_linkXYZ = size(ratio4_linkXYZ(ratio4_linkXYZ<=threshold1(1)),1)/size(color4_linkXYZ,1);
cross42_linkXYZ=size(ratio4_linkXYZ(ratio4_linkXYZ>threshold1(2)&ratio4_linkXYZ<=threshold2(1)),1)/size(color4_linkXYZ,1);
cross43_linkXYZ=size(ratio4_linkXYZ(ratio4_linkXYZ>threshold2(2)&ratio4_linkXYZ<=threshold3(1)),1)/size(color4_linkXYZ,1);
cross44_linkXYZ=size(ratio4_linkXYZ(ratio4_linkXYZ>threshold3(2)),1)/size(color4_linkXYZ,1);


goodness1_linkXYZ = size(color1_linkXYZ,1)/size(ratio1_linkXYZ,1);
goodness2_linkXYZ = size(color2_linkXYZ,1)/size(ratio2_linkXYZ,1);
goodness3_linkXYZ = size(color3_linkXYZ,1)/size(ratio3_linkXYZ,1);
goodness4_linkXYZ = size(color4_linkXYZ,1)/size(ratio4_linkXYZ,1);
goodnessALL_linkXYZ =[goodness1_linkXYZ;goodness2_linkXYZ;goodness3_linkXYZ;goodness4_linkXYZ];




% Individual Fit**************************************************************************************************************************************************************
[P1up,CRLB1up, LL1up] =  mleFit_LM(output1up,5,50,coeff1,0,1,zstart(1,:));
[P1down,CRLB1down, LL1down] =  mleFit_LM(output1down,5,50,coeff2,0,1,zstart(1,:));

[P2up,CRLB2up, LL2up] =  mleFit_LM(output2up,5,50,coeff1,0,1,zstart(1,:));
[P2down,CRLB2down, LL2down] =  mleFit_LM(output2down,5,50,coeff2,0,1,zstart(1,:));

[P3up,CRLB3up, LL3up] =  mleFit_LM(output3up,5,50,coeff1,0,1,zstart(1,:));
[P3down,CRLB3down, LL3down] =  mleFit_LM(output3down,5,50,coeff2,0,1,zstart(1,:));

[P4up,CRLB4up, LL4up] =  mleFit_LM(output4up,5,50,coeff1,0,1,zstart(1,:));
[P4down,CRLB4down, LL4down] =  mleFit_LM(output4down,5,50,coeff2,0,1,zstart(1,:));

% calculate ratio
ratio1_sep = (P1up(:,3)-P1down(:,3))./(P1up(:,3)+P1down(:,3));
ratio2_sep = (P2up(:,3)-P2down(:,3))./(P2up(:,3)+P2down(:,3));
ratio3_sep = (P3up(:,3)-P3down(:,3))./(P3up(:,3)+P3down(:,3));
ratio4_sep = (P4up(:,3)-P4down(:,3))./(P4up(:,3)+P4down(:,3));
ratio_sep_ALL = [ratio1_sep;ratio2_sep;ratio3_sep;ratio4_sep];

%define threshold
threshold1 = 0.5457;
threshold2 = 0.7668;
threshold3 = 0.9116;

IntervalT = 0.01;
threshold1 = [0.5457-IntervalT 0.5457+IntervalT];
threshold2 = [0.7668-IntervalT 0.7668+IntervalT];
threshold3 = [0.9116-IntervalT 0.9116+IntervalT];

% separaration of individualfit
color1_separate = ratio1_sep(ratio1_sep<=threshold1(1)|(ratio1_sep>threshold1(2)&ratio1_sep<=threshold2(1))|(ratio1_sep>threshold2(2)&ratio1_sep<=threshold3(1)|ratio1_sep>threshold3(2)));
cross11_separate = size(ratio1_sep(ratio1_sep<=threshold1(1)),1)/size(color1_separate,1);
cross12_separate=size(ratio1_sep(ratio1_sep>threshold1(2)&ratio1_sep<=threshold2(1)),1)/size(color1_separate,1);
cross13_separate=size(ratio1_sep(ratio1_sep>threshold2(2)&ratio1_sep<=threshold3(1)),1)/size(color1_separate,1);
cross14_separate=size(ratio1_sep(ratio1_sep>threshold3(2)),1)/size(color1_separate,1);

color2_separate = ratio2_sep(ratio2_sep<=threshold1(1)|(ratio2_sep>threshold1(2)&ratio2_sep<=threshold2(1))|(ratio2_sep>threshold2(2)&ratio2_sep<=threshold3(1)|ratio2_sep>threshold3(2)));
cross21_separate = size(ratio2_sep(ratio2_sep<=threshold1(1)),1)/size(color2_separate,1);
cross22_separate=size(ratio2_sep(ratio2_sep>threshold1(2)&ratio2_sep<=threshold2(1)),1)/size(color2_separate,1);
cross23_separate=size(ratio2_sep(ratio2_sep>threshold2(2)&ratio2_sep<=threshold3(1)),1)/size(color2_separate,1);
cross24_separate=size(ratio2_sep(ratio2_sep>threshold3(2)),1)/size(color2_separate,1);

color3_separate = ratio3_sep(ratio3_sep<=threshold1(1)|(ratio3_sep>threshold1(2)&ratio3_sep<=threshold2(1))|(ratio3_sep>threshold2(2)&ratio3_sep<=threshold3(1)|ratio3_sep>threshold3(2)));
cross31_separate = size(ratio3_sep(ratio3_sep<=threshold1(1)),1)/size(color3_separate,1);
cross32_separate=size(ratio3_sep(ratio3_sep>threshold1(2)&ratio3_sep<=threshold2(1)),1)/size(color3_separate,1);
cross33_separate=size(ratio3_sep(ratio3_sep>threshold2(2)&ratio3_sep<=threshold3(1)),1)/size(color3_separate,1);
cross34_separate=size(ratio3_sep(ratio3_sep>threshold3(2)),1)/size(color3_separate,1);

color4_separate = ratio4_sep(ratio4_sep<=threshold1(1)|(ratio4_sep>threshold1(2)&ratio4_sep<=threshold2(1))|(ratio4_sep>threshold2(2)&ratio4_sep<=threshold3(1)|ratio4_sep>threshold3(2)));
cross41_separate = size(ratio4_sep(ratio4_sep<=threshold1(1)),1)/size(color4_separate,1);
cross42_separate=size(ratio4_sep(ratio4_sep>threshold1(2)&ratio4_sep<=threshold2(1)),1)/size(color4_separate,1);
cross43_separate=size(ratio4_sep(ratio4_sep>threshold2(2)&ratio4_sep<=threshold3(1)),1)/size(color4_separate,1);
cross44_separate=size(ratio4_sep(ratio4_sep>threshold3(2)),1)/size(color4_separate,1);


goodness1_sep = size(color1_separate,1)/size(ratio1_sep,1);
goodness2_sep = size(color2_separate,1)/size(ratio2_sep,1);
goodness3_sep = size(color3_separate,1)/size(ratio3_sep,1);
goodness4_sep = size(color4_separate,1)/size(ratio4_sep,1);
goodnessALL_sep =[goodness1_sep;goodness2_sep;goodness3_sep;goodness4_sep];


%% export
% crosstalk of the 4 dyes under different fitting schemes
cross_linkxyzP = [cross11_linkPF cross12_linkPF cross13_linkPF cross14_linkPF;cross21_linkPF cross22_linkPF cross23_linkPF cross24_linkPF;cross31_linkPF cross32_linkPF cross33_linkPF cross34_linkPF;cross41_linkPF cross42_linkPF cross43_linkPF cross44_linkPF];
cross_linkxyz = [cross11_linkXYZ cross12_linkXYZ cross13_linkXYZ cross14_linkXYZ;cross21_linkXYZ cross22_linkXYZ cross23_linkXYZ cross24_linkXYZ;cross31_linkXYZ cross32_linkXYZ cross33_linkXYZ cross34_linkXYZ;cross41_linkXYZ cross42_linkXYZ cross43_linkXYZ cross44_linkXYZ];
cross_separate = [cross11_separate cross12_separate cross13_separate cross14_separate;cross21_separate cross22_separate cross23_separate cross24_separate;cross31_separate cross32_separate cross33_separate cross34_separate;cross41_separate cross42_separate cross43_separate cross44_separate];



