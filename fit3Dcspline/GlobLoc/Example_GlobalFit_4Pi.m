% Copyright (c) 2021 Li Lab, Southern University of Science and Technology, Shenzhen & Ries Lab, European Molecular Biology Laboratory, Heidelberg.
% author: Sheng Liu & Yiming Li 
% email: liym2019@sustech.edu.cn & sheng.liu@embl.de
% date: 2021.08.27
% Tested with CUDA 11.3 (Express installation) and Matlab 2019a
%%
clearvars;
resultfolder = 'output\';
if ~exist(resultfolder,'dir')
    mkdir(resultfolder)
end

%% prepare parameters
paraSim.NA = 1.35;                                                % numerical aperture of obj             
paraSim.refmed = 1.406;                                            % refractive index of sample medium
paraSim.refcov = 1.518;                                           % refractive index of converslip
paraSim.refimm = 1.406;                                           % refractive index of immersion oil
paraSim.lambda = 668;                                             % wavelength of emission
paraSim.objStage0_upper = -0;                                        % nm, initial objStage0 position,relative to focus at coverslip
paraSim.objStage0_lower = -0;                                        % nm, initial objStage0 position,relative to focus at coverslip
paraSim.zemit0_upper = -1*paraSim.refmed/paraSim.refimm*(paraSim.objStage0_upper);  % reference emitter z position, nm, distance of molecule to coverslip
paraSim.zemit0_lower = -1*paraSim.refmed/paraSim.refimm*(paraSim.objStage0_lower);  % reference emitter z position, nm, distance of molecule to coverslip
paraSim. pixelSizeX = 120;                                        % nm, pixel size of the image
paraSim. pixelSizeY = 120;                                        % nm, pixel size of the image
paraSim.Npupil = 64;                                             % sampling at the pupil plane
paraSim.aberrations(:,:,1) = [2,-2,0.0; 2,2,-0.06; 3,-1,0.0; 3,1,0.0; 4,0,0; 3,-3,0.0; 3,3,0.0; 4,-2,0.0; 4,2,0.00; 5,-1,0.0; 5,1,0.0; 6,0,0.0; 4,-4,0.0; 4,4,0.0;  5,-3,0.0; 5,3,0.0;  6,-2,0.0; 6,2,0.0; 7,1,0.0; 7,-1,0.00; 8,0,0.0];
paraSim.aberrations(:,:,2) = [2,-2,0.0; 2,2,0.06; 3,-1,0.0; 3,1,0.0; 4,0,0.0; 3,-3,0.0; 3,3,0.0; 4,-2,0.0; 4,2,0.00; 5,-1,0.0; 5,1,0.0; 6,0,0.0; 4,-4,0.0; 4,4,0.0;  5,-3,0.0; 5,3,0.0;  6,-2,0.0; 6,2,0.0; 7,1,0.0; 7,-1,0.00; 8,0,0.0];
paraSim.aberrations(:,3,:) =  paraSim.aberrations(:,3,:)*paraSim.lambda;
paraSim.offset = [0 0];
paraSim.phaseshift = [0 ,pi/2, pi, 3*pi/2];

Nmol = 101;
Npixels = 25;
paraSim.Nmol = Nmol;
paraSim.sizeX = Npixels;
paraSim.sizeY = Npixels;

paraSim.xemit = zeros(1,Nmol);                             %nm
paraSim.yemit = zeros(1,Nmol);                             %nm
paraSim.zemit = linspace(-1000,1000,Nmol);                 %nm
paraSim.objStage = zeros(1,Nmol);                          %nm
%% simulate IAB model
% simulate vectorial PSFs
[PSFs,PSFsUpper,PSFsLower,WaberrationUpper, WaberrationLower] = vectorPSF_4Pi(paraSim);

% calculate IAB model
ipalm_im  = PSFs;
phaseshift = paraSim.phaseshift;
k = 2 * pi / paraSim.lambda; %lambda in nm
zcand = paraSim.zemit;% if move obj stage use paraSim.objStage
zstep = zcand(2) - zcand(1);
imsz = paraSim.sizeX;
I = squeeze((ipalm_im(:, :, :, 1) + ipalm_im(:, :, :, 3)) / 2);
kz2 = 2 * k * zcand';
kz2 = permute(repmat(kz2, 1, imsz, imsz), [2, 3, 1]);
F_phi1 = squeeze(ipalm_im(:, :, :, 1)) - I;
F_phi2 = squeeze(ipalm_im(:, :, :, 2)) - I;
phi1 = phaseshift(1);
phi2 = phaseshift(2);
A = (F_phi1 .* sin(kz2 + phi2) - F_phi2 .* sin(kz2 + phi1)) / sin(phi2 - phi1);
B = (-F_phi1 .* cos(kz2 + phi2) + F_phi2 .* cos(kz2 + phi1)) / sin(phi2 - phi1);
phi4 = phaseshift(4);
check_PSF = I + A .* cos(kz2 + phi4) + B .* sin(kz2 + phi4);
Ispline = Spline3D_interp(I./2);
Aspline = Spline3D_interp(A./2);
Bspline = Spline3D_interp(B./2);
IAB = cat(5,Aspline,Bspline,Ispline);
IABall = cat(6,IAB,IAB,IAB,IAB); % assume all quadrants have the same IAB with only a constant phase shift, which will be added for data simulation and localization 

%% simulate data
zT = paraSim.lambda/zstep/2;
bxsz = 15;
Nz = 121;
Nrep = 100;
z0 = repmat(linspace(26,76,Nz)',Nrep,1);
Nfit = Nz*Nrep;
x0 = bxsz/2-0.5+rand(Nfit,1)*2-1;
y0 = bxsz/2-0.5+rand(Nfit,1)*2-1;
phi0 = z0./zT.*2*pi;
I0 = ones(Nfit,1);
bg0 = zeros(Nfit,1);
psfsim = zeros(bxsz,bxsz,Nfit,4);
sharedA = ones(6,Nfit);
phi0A = ones(4,Nfit).*[0 ,pi/2, pi, 3*pi/2]';
dTAll = zeros(6,4,Nfit);
iterations = 0;
[~,~,~,psfh] = GPUmleFit_LM_4Pi_v1(single(psfsim),uint32(sharedA),iterations,single(IABall),single(dTAll),single(phi0A),single(z0),single(phi0),single(y0),single(x0),single(I0),single(bg0));
psfsim = reshape(psfh,bxsz,bxsz,[],4);
photon = 2000;
bg = 20;
psf = psfsim*photon + bg;
data = single(poissrnd(psf)); 

%% calculate CRLB
% unlink
sharedA = repmat([1,1,0,0,1,1]',[1 Nfit]);
[~,ch_uk] = GPUmleFit_LM_4Pi_v1(single(psfsim),uint32(sharedA),iterations,single(IABall),single(dTAll),single(phi0A),single(z0),single(phi0),single(y0),single(x0),single(I0.*photon),single(bg0+bg));
stdM_uk = reshape(sqrt(ch_uk),[],1,12);
% link
sharedA = repmat([1,1,1,1,1,1]',[1 Nfit]);
[ph,ch, Lh, psfh] = GPUmleFit_LM_4Pi_v1(single(psfsim),uint32(sharedA),iterations,single(IABall),single(dTAll),single(phi0A),single(z0),single(phi0),single(y0),single(x0),single(I0.*photon),single(bg0+bg));
stdM = reshape(sqrt(ch),[],1,6);

%% localization
% unlink
N = size(data,3);
sr.BGoffset = 0;
sr.Iratio = [1,1,1,1];
sr.Phi0 = [0 ,pi/2, pi, 3*pi/2];
sr.Boxsize = bxsz;
sr.Initz = size(A,3)/2;
sr.InitPhase = [0,pi/3,2*pi/3];
sr.Initx = [];
sr.Dz = [0,0,0,0];
sr.Dphi = [0,0,0,0];
sr.IABall = IABall;
shared = [1,1,0,0,1,1];
dx = zeros(N,4);
dy = zeros(N,4);
[P_uk] = psfloc(sr,data,dx,dy,shared);

% link
N = size(data,3);
sr.BGoffset = 0;
sr.Iratio = [1,1,1,1];
sr.Phi0 = [0 ,pi/2, pi, 3*pi/2];
sr.Boxsize = bxsz;
sr.Initz = size(A,3)/2;
sr.InitPhase = [0,pi/3, 2*pi/3];
sr.Initx = [];
sr.Dz = [0,0,0,0];
sr.Dphi = [0,0,0,0];
sr.IABall = IABall;
shared = [1,1,1,1,1,1];
dx = zeros(N,4);
dy = zeros(N,4);
[P] = psfloc(sr,data,dx,dy,shared);

% photometry
pz = zstep;
datasum = squeeze(sum(data,4));
datasum(datasum<=0) = 1e-6;
PSFsigma = 1.15;
iterations = 100;
fittype = 4;
[P_pm] = mleFit_LM(single(datasum),fittype,iterations,[PSFsigma,PSFsigma]); 
xf = P_pm(:,2);
yf = P_pm(:,1);
[xx,yy] = meshgrid(0:bxsz-1,0:bxsz-1);
comx = squeeze(sum(datasum.*xx,[1,2]))./squeeze(sum(datasum,[1,2]));
comy = squeeze(sum(datasum.*xx,[1,2]))./squeeze(sum(datasum,[1,2]));
maskx = xf>bxsz-5|xf<4;
masky = yf>bxsz-5|yf<4;
xf(maskx)=comx(maskx);
yf(masky)=comy(masky);
zT_ast = zT*pz;
phi_s = 0;
phi_p = pi/2;
[z_ang,z_phi_a,z_phi,phi] = phifitM(data,xf,yf,phi_s,phi_p,zT_ast,0);

%% calculate localization precision and bias
% unlink
pz = zstep;
pxsz = paraSim.pixelSizeX;
dphi = wrapToPi(unwrap(reshape(P_uk(:,12)-z0./zT.*2*pi,Nz,[])));
dphi = dphi(:).*zT*pz/2/pi;
dz = (reshape(P_uk(:,11)-z0,Nz,[])).*pz;
dz = dz(:);
dx = (reshape(P_uk(:,2)-x0,Nz,[])).*pxsz;
dy = (reshape(P_uk(:,1)-y0,Nz,[])).*pxsz;
dx = dx(:);
dy = dy(:);
xyzp = cat(3,dx,dy,dphi,dz);
mask = abs(dz)>50 | abs(dx)>100 | abs(dy)>100 | abs(dphi)>50;
xyzp(mask,:) = 0;
xyzp_mean_uk = squeeze(sum(reshape(xyzp,Nz,[],4),2))./squeeze(sum(reshape(~mask,Nz,[],1),2));
xyzp_rms_uk = squeeze(sqrt(sum(reshape(xyzp,Nz,[],4).^2,2)./sum(reshape(~mask,Nz,[],1),2)));
xyzp_crlb_uk = squeeze(mean(reshape(stdM_uk(:,:,[2,1,12,11]),Nz,[],4),2)).*[pxsz,pxsz,zT*pz/2/pi,pz];
xyzp_std_uk = squeeze(std(reshape(xyzp,Nz,[],4),1,2));

% link
pz = zstep;
pxsz = paraSim.pixelSizeX;
dphi = wrapToPi(unwrap(reshape(P(:,6)-z0./zT.*2*pi,Nz,[])));
dphi = dphi(:).*zT*pz/2/pi;
dz = (reshape(P(:,5)-z0,Nz,[])).*pz;
dz = dz(:);
dx = (reshape(P(:,2)-x0,Nz,[])).*pxsz;
dy = (reshape(P(:,1)-y0,Nz,[])).*pxsz;
dx = dx(:);
dy = dy(:);
xyzp = cat(3,dx,dy,dphi,dz);
mask = abs(dz)>100 | abs(dx)>100 | abs(dy)>100;
xyzp(mask,:) = 0;
xyzp_mean = squeeze(sum(reshape(xyzp,Nz,[],4),2))./squeeze(sum(reshape(~mask,Nz,[],1),2));
xyzp_rms = squeeze(sqrt(sum(reshape(xyzp,Nz,[],4).^2,2)./sum(reshape(~mask,Nz,[],1),2)));
xyzp_crlb = squeeze(mean(reshape(stdM(:,:,[2,1,6,5]),Nz,[],4),2)).*[pxsz,pxsz,zT*pz/2/pi,pz];
xyzp_std = squeeze(std(reshape(xyzp,Nz,[],4),1,2));

% photometry
dphi = wrapToPi(reshape(phi-(z0-Nmol/2)./zT.*2*pi,Nz,[]));
dphi = (dphi(:)-median(dphi(:))).*zT_ast/2/pi;
dx = (reshape(xf-x0,Nz,[])).*pxsz;
dy = (reshape(yf-y0,Nz,[])).*pxsz;
dx = dx(:);
dy = dy(:);
xyp = cat(3,dx,dy,dphi);
mask = abs(xyp)>250;
xyp(mask) = 0;
xyp_mean_pm = squeeze(sum(reshape(xyp,Nz,[],3),2))./squeeze(sum(reshape(~mask,Nz,[],3),2));
xyp_rms_pm = squeeze(sqrt(sum(reshape(xyp,Nz,[],3).^2,2)./sum(reshape(~mask,Nz,[],3),2)));
xyp_std_pm = squeeze(std(reshape(xyp,Nz,[],3),1,2));

%% plot result
zs = linspace(-500,500,Nz)';
label = {'x','y','z'};
h = figure('Position',[200,300,1310,270]);

for ii = 1:3
    ha = subplot(1,4,ii);hold on
    plot(zs,xyzp_crlb(:,ii),'k')
    plot(zs(1:4:end),xyp_std_pm(1:4:end,ii),'^','markersize',3)
    plot(zs(1:4:end),xyzp_std(1:4:end,ii),'o','markersize',3)
    ha.Box = 'on';
    ha.XLabel.String = 'z (nm)';
    ha.YLabel.String = [label{ii}, ' loc prec (nm)'];
end
plot(zs,xyzp_crlb_uk(:,ii),'color',[1,1,1].*0.5)
legend('CRLB','photometry','Model Fit','CRLB unlink')

ha = subplot(1,4,4);hold on
plot(zs,zeros(size(zs)),'k')
plot(zs(1:4:end),xyp_mean_pm(1:4:end,3),'^','markersize',3)
plot(zs(1:4:end),xyzp_mean(1:4:end,3),'o','markersize',3)
ha.XLabel.String = 'z (nm)';
ha.YLabel.String = 'bias in z (nm)';
ha.Box = 'on';
%% save figure
set(h,'PaperPositionMode','auto')
imgname = 'std_crlb_vec';
print(h,'-dpng','-r300',[resultfolder,imgname])
print(h,'-depsc','-r300',[resultfolder,imgname])




