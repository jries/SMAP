%  Copyright (c)2017 Ries Lab, European Molecular Biology Laboratory,
%  Heidelberg.
%  
%  This file is part of GPUmleFit_LM Fitter.
%  
%  GPUmleFit_LM Fitter is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%  
%  GPUmleFit_LM Fitter is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%  
%  You should have received a copy of the GNU General Public License
%  along with GPUmleFit_LM Fitter.  If not, see <http://www.gnu.org/licenses/>.
%  
%  
%  Additional permission under GNU GPL version 3 section 7
%  
%  If you modify this Program, or any covered work, by
%  linking or combining it with libraries required for interaction
%  with analysis programs such as Igor Pro or Matlab,
%  the licensors of this Program grant you additional permission
%  to convey the resulting work.

%% example script for fitting of 3D data
%Installation instructions
%  Requires Matlab 2016a or newer including the following tool boxes:
%  image analysis, curve fitting, optimization, statistics
%  GPU support for CUDA 8


%% add path to helper functions
% addpath('bfmatlab') %now replaced by own tiff reader mytiffread.m
addpath('shared')

%% make bead calibration
%run 3D calibration GUI (alternatively, you can directly call calibrate3D with proper parameters)
%and make 3D calibration
%e.g. using bead files: in data/beadstacks_3D_astig (extracted beadstacks_3D_astig.zip)
% save e.g. as data/bead_3dcal.mat
%For fitting 2D datasets, create 3D calibration from 2D PSF stack (e.g. data/beadstacks_2D)
% save e.g. as data/bead2D_3dcal.mat

g=calibrate3D_GUI;  %look at the tooltips for explanantions of the parameters
fill=[fileparts(pwd) filesep 'data/bead_3dcal.mat']; %refer to data outside the main directory
g.guihandles.outputfile.String=fill;
g.guihandles.filelist.String=[fileparts(pwd) filesep 'data/beadstacks_3D_astig/stack3D_1.tif']; %you can load additional bead stacks in the GUI
%% load bead calibration
cal=load(fill); %load bead calibration


%% either simulate data or load experimental tiff file
mode = 1;
if mode ==1 % simulate data
    p.offset=0;
    p.conversion=1;
    
    %numlocs: number of simulated molecules. For maximum fitting speeds on  the GPU, 
    %it should be >10^5. For smaller GPUs, adust to avoid memory overflow
    %error to 10^5-10^5. For CPU fitter, use 10^3-10^4.
    numlocs=50000;
    RoiPixelsize=13; %ROI in pixels for simulation
    dz=cal.cspline.dz;  %coordinate system of spline PSF is corner based and in units pixels / planes
    z0=cal.cspline.z0; % distance and midpoint of stack in spline PSF, needed to translate into nm coordinates
    dx=floor(RoiPixelsize/2); %distance of center from edge
    ground_truth.z=linspace(-500,500,numlocs)'; %define some coordinates. Alternatively, use rand
    ground_truth.x=linspace(-0.5,0.5,numlocs)';
    ground_truth.y=sin(ground_truth.x*4*pi);
    coordinates=horzcat(ground_truth.x+dx,ground_truth.y+dx,ground_truth.z/dz+z0);
    Intensity=2000; %number of photons / localization
    background=10; %number of background photons per pixel
    imstack = simSplinePSF(RoiPixelsize,cal.cspline.coeff,Intensity,background,coordinates); %simulate images
    
else %experimental data
    
    file=[fileparts(pwd) filesep 'data/single_bead.tif']; %experimental astgmatic bead stack.
    imstackadu=readfile_tif(file); % Stack of ROIs in photons.
    p.offset=400;p.conversion=.1;
    imstack=(single(imstackadu)-p.offset)/p.conversion;% if necessary, convert ADU into photons.
    ground_truth.z=((1:size(imstack,3))'-size(imstack,3)/2+1)*10; %dz=10 nm; convert frame to nm   
end
%% load sCMOS varmap or set sCMOSvarmap to zero
% sCMOSvarmap=ones(size(imstack),'single'); %the variance map, same size as imstack, single format. 
sCMOSvarmap=0; %if scalar: use EMCCD fitter;

%% fit image stacks
numlocs=size(imstack,3); 
imstack=single(imstack); %imstack needs to be in photons. The fitters require the stacks in single-format;

if isfield(cal,'gauss_zfit')&&~isempty(cal.gauss_zfit) %only if Gaussian PSF model was calibrated in the calibrate3D_GUI, we can test the Gaussian fit
    %fit Gaussian model, direct z fit, emCCD mode
    tic
    [P,CRLB]=mleFit_LM(imstack,3,50,single(cal.gauss_zfit),sCMOSvarmap,1);
    tgz=toc;
    disp(['Gauss z: ' num2str(numlocs/tgz) ' fits/s']);
    gaussz.x=P(:,1);gaussz.y=P(:,2); %x,y in pixels 
    gaussz.z=P(:,5)*1000;


    %fit elliptical Gaussian model, extract z from sx, sy. emCCD mode
    tic
    [P,CRLB]=mleFit_LM(imstack,4,50,1,sCMOSvarmap,1);
    tgsxsy=toc;
    disp(['Gaussxy: ' num2str(numlocs/tgsxsy) ' fits/s']);

    gausssxsy.x=P(:,1);gausssxsy.y=P(:,2); %x,y in pixels 
    sx=P(:,5);sy=P(:,6);
    gausssxsy.z=sxsy2z(sx,sy,cal.gauss_sx2_sy2); %z from sx, sy

else
    gaussz.x=zeros(size(imstack,3),1);gaussz.y=zeros(size(imstack,3),1); gaussz.z=zeros(size(imstack,3),1);
    gausssxsy.x=zeros(size(imstack,3),1);gausssxsy.y=zeros(size(imstack,3),1); gausssxsy.z=zeros(size(imstack,3),1);
end

%fit experimental astigmatic PSF, cspline, emCCD mode
tic
[Pcspline,CRLB]=mleFit_LM(imstack,5,50,single(cal.cspline.coeff),sCMOSvarmap,1);
tspline=toc;
disp(['cspline: ' num2str(numlocs/tspline) ' fits/s']);
dx=floor(size(imstack,1)/2);
cspline.x=Pcspline(:,1)-dx;cspline.y=Pcspline(:,2)-dx; %x,y in pixels 
cspline.z=(Pcspline(:,5)-cal.cspline.z0)*cal.cspline.dz;

z_err=sqrt(CRLB(:,5))*cal.cspline.dz; %CRLB is the variance of the parameters, to obtain a localization precision we have to take the square root

% calculate error for all fits. 
zrange=400; %Only take into accoutn central part. Focus +/- zrange
inz=abs(ground_truth.z)<zrange;

%difference between fitted z and ground truth
cspline.zrel=cspline.z-ground_truth.z;
gaussz.zrel=gaussz.z-ground_truth.z;
gaussz.zrel(isnan(gaussz.zrel))=0;
gausssxsy.zrel=gausssxsy.z-ground_truth.z;
gausssxsy.zrel(isnan(gausssxsy.zrel))=0;

%numpoints
numpoints=min(2000,size(imstack,3));
range=1:round(length(ground_truth.z)/numpoints):length(ground_truth.z);
%plot fitted z vs ground-truth z
    figure(101)
    hold off
    plot(ground_truth.z(range),gaussz.zrel(range)-mean(gaussz.zrel(inz)),'r.')
    hold on
    plot(ground_truth.z(range),gausssxsy.zrel(range)-mean(gausssxsy.zrel(inz)),'.','Color',[0.7,0.7,0])
    plot(ground_truth.z(range),cspline.zrel(range)-mean(cspline.zrel(inz)),'b.')
    gs=fit(ground_truth.z(range),double(gaussz.zrel(range)-mean(gaussz.zrel(inz))),'smoothingspline','SmoothingParam',.00001);
    gss=fit(ground_truth.z(range),double(gausssxsy.zrel(range)-mean(gausssxsy.zrel(inz))),'smoothingspline','SmoothingParam',.00001);
    gz=fit(ground_truth.z(range),double(cspline.zrel(range)-mean(cspline.zrel(inz))),'smoothingspline','SmoothingParam',.00001);

    plot(ground_truth.z(range),ground_truth.z(range)*0,'k','LineWidth',3)
    plot(ground_truth.z(range),gs(ground_truth.z(range)),'r','LineWidth',3)
    plot(ground_truth.z(range),gss(ground_truth.z(range)),'Color',[0.7,0.7,0],'LineWidth',3)
    plot(ground_truth.z(range),gz(ground_truth.z(range)),'b','LineWidth',3)

    plot([-zrange -zrange],[-1000 1000],'c')
    plot([zrange zrange],[-1000 1000],'c')

    ss=std(cspline.zrel(inz))*10;
    ylim([-ss ss])
    xlabel('ground truth z (nm)')
    ylabel('z (fitted) - z (ground truth) (nm)')

    title('Accuracy')

    cspline.dz=std(ground_truth.z(inz)-cspline.z(inz),'omitnan');
    zgauss.dz=std(ground_truth.z(inz)-gaussz.z(inz),'omitnan');
    gausssxsy.dz=std(ground_truth.z(inz)-gausssxsy.z(inz),'omitnan');

    legendtxt{4}='smoothing spline';
    legendtxt{3}=['spline fit. error: ' num2str(cspline.dz, '%3.1f') ' nm'];
    legendtxt{1}=['Gaussian z fit. error: ' num2str(zgauss.dz, '%3.1f') ' nm'];
    legendtxt{2}=['Gaussian sx, sy fit. error: ' num2str(gausssxsy.dz, '%3.1f') ' nm'];
    legend(legendtxt)

%calculate lateral error
if mode == 1
    cspline.dx=std(ground_truth.x(inz)-cspline.x(inz),'omitnan'); %in pixels
    cspline.dy=std(ground_truth.y(inz)-cspline.y(inz),'omitnan');
    dx_gaussz=std(ground_truth.x(inz)-gaussz.x(inz),'omitnan');
    dy_gaussz=std(ground_truth.y(inz)-gaussz.y(inz),'omitnan');
    gausssxsy.dx=std(ground_truth.x(inz)-gausssxsy.x(inz),'omitnan');
    gausssxsy.dy=std(ground_truth.y(inz)-gausssxsy.y(inz),'omitnan');
    %plot 3D scatter plot
   
    figure(102);
    hold off
    scatter3(ground_truth.x(range),ground_truth.y(range),ground_truth.z(range),3);
    hold on
    scatter3(cspline.x(range),cspline.y(range),cspline.z(range),5);
    xlabel('x (nm)');ylabel('y (nm)');zlabel('z (nm)');
    legend('ground truth','cspline fit')
    title('fitted vs ground truth positions')
else
    cspline.dx=std(cspline.x(inz),'omitnan'); %in pixels
    cspline.dy=std(cspline.y(inz),'omitnan');
    dx_gaussz=std(gaussz.x(inz),'omitnan');
    dy_gaussz=std(gaussz.y(inz),'omitnan');
    gausssxsy.dx=std(gausssxsy.x(inz),'omitnan');
    gausssxsy.dy=std(gausssxsy.y(inz),'omitnan');
end
%% make 2D PSF
g=calibrate3D_GUI;  %look at the tooltips for explanantions of the parameters
fill=[fileparts(pwd) filesep 'data/bead2d_3dcal.mat']; %refer to data outside the main directory
g.guihandles.outputfile.String=fill;
g.guihandles.filelist.String=[fileparts(pwd) filesep 'data/beadstacks_2D/stack2D_1.tif'];
%% fit 2D dataset with cspline
% Here we simulate Nfits positions per data point and calculate the
% z-dependent error
cal=load([fileparts(pwd) filesep 'data/bead2d_3dcal.mat']); %load bead calibration

ztruth = -475:50:475; %z positions for which we want to simulate fluorophores
Nfits = 300; %fits per data points
Nphotons = 2000; %photons/localizations
Npixels = 17; %size of the ROI
bg = 10; %bg photons per pixel
dx = 132; %nm pixel size
dy = 132; %nm
dz=cal.cspline.dz;  %coordinate system of spline PSF is corner based and in units pixels / planes
z0=cal.cspline.z0; % distance and midpoint of stack in spline PSF, needed to translate into nm coordinates
clear Cspline2DF fractionmisassigned
for i = 1: length(ztruth)
    disp([num2str(i) ' of ' num2str(length(ztruth))])
    coordsxy = Npixels/2 -1 +2*rand([Nfits 2]); %random x, y positions
    coordsz = ztruth(i)/dz+z0*ones(Nfits,1);
    coordinates = [coordsxy coordsz];
    imstack = simSplinePSF(Npixels,cal.cspline.coeff,Nphotons,bg,coordinates);
    zstartnm=400; %z start parameter 400 nm above and below the focus
    zstart=[-1 1]*zstartnm/cal.cspline.dz; %in units of dz
    
    [P CRLB LL ]=mleFit_LM(imstack,5,50,single(cal.cspline.coeff),sCMOSvarmap,1,zstart); 
    
    
    z=(P(:,5)-z0).*dz;
    locprecz=sqrt(CRLB(:,5))*dz;
    %determine fraction of misassigned localizations
    misassigned=sign(ztruth(i))~=sign(z) & abs(z)> locprecz; %those close to the focus cannot be distinguished and scatter around zero. By excluding only those within their localization precision, we probably overestimate the misassignments.
    fractionmisassigned(i)=sum(misassigned)/length(misassigned);
    
    %Filter
    %to calculate localization accuracy and precision we here exclude the
    %misassigned localizations, as they result in very large effective
    %errors, although in practice they form a mirror image of the structure
    %and can usually be neglected.
    Cspline2DF.x = P(~misassigned,1);
    Cspline2DF.y = P(~misassigned,2);
    Cspline2DF.z= (P(~misassigned,5)-z0).*dz;
    
    Cspline2DF.CRLBx = CRLB(~misassigned,1);
    Cspline2DF.CRLBy = CRLB(~misassigned,2);
    Cspline2DF.CRLBz = CRLB(~misassigned,5);
    
    Cspline2DF.s_x_found(i,1) = std( Cspline2DF.x-coordsxy(~misassigned,1));
    Cspline2DF.s_y_found(i,1) = std( Cspline2DF.y-coordsxy(~misassigned,2));
    Cspline2DF.s_z_found(i,1) = std( Cspline2DF.z-ztruth(i));
    
    Cspline2DF.meansqrtCRLBx(i,1) = mean(sqrt(Cspline2DF.CRLBx));
    Cspline2DF.meansqrtCRLBy(i,1) = mean(sqrt(Cspline2DF.CRLBy));
    Cspline2DF.meansqrtCRLBz(i,1) = mean(sqrt(Cspline2DF.CRLBz))*dz;
    
    
    Cspline2DF.RMSEX(i,1) = sqrt(mean((Cspline2DF.x-coordsxy(~misassigned,1)).^2));
    Cspline2DF.RMSEY(i,1) = sqrt(mean((Cspline2DF.y-coordsxy(~misassigned,2)).^2));
    Cspline2DF.RMSEZ(i,1) = sqrt(mean((Cspline2DF.z-ztruth(i)).^2));

end

%Plot localization precision and fraction of misassigned localizations
    figure(103);subplot(1,2,1)
    hold off
    plot(ztruth,Cspline2DF.meansqrtCRLBx*dx,'-');
    hold on
    plot(ztruth,Cspline2DF.meansqrtCRLBy*dy,'-');
    plot(ztruth,Cspline2DF.meansqrtCRLBz,'-');
    plot(ztruth,Cspline2DF.s_x_found*dx,'o');
    plot(ztruth,Cspline2DF.s_y_found*dy,'o');
    plot(ztruth,Cspline2DF.s_z_found,'o');
    xlabel('z (nm)')
    ylabel('x,y,z localization precision (nm)')
    legend('CRLBx','CRLBy','CRLBz','std(x)','std(y)','std(z)');
    title('localization precision 2D PSF');
    subplot(1,2,2)
    plot(ztruth, fractionmisassigned,'-o')
    title('fraction of misassigned localizations')
    xlabel('z (nm)')
    ylabel('fraction')
