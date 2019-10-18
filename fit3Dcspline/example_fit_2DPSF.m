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

%% example script for fitting of 2D data
%Installation instructions
%  Requires Matlab 2016a or newer including the following tool boxes:
%  image analysis, curve fitting, optimization, statistics
%  GPU support for CUDA 8


%% add path to helper functions
addpath('shared')

%% make bead calibration
%run 3D calibration GUI and make 3D calibration
%For fitting 2D datasets, create 3D calibration from 2D PSF stack (e.g. example_data/beadstacks_2D)
% save e.g. as data/bead_2DPSF_3dcal.mat

% or generate calibration files programmatically:
if ~exist([pwd filesep 'example_data' filesep 'bead_2DPSF_3dcal.mat'],'file') %only run if no calibration files have been generated, as this takes time.
    example_calibration
end


%% load bead calibration
cal=load([pwd filesep 'example_data' filesep 'bead_2DPSF_3dcal.mat']); %load bead calibration

%% fit 2D dataset with cspline
% Here we simulate Nfits positions per data point and calculate the
% z-dependent error

ztruth = -725:50:625; %z positions for which we want to simulate fluorophores
Nfits = 300; %fits per data points
Nphotons = 3000; %photons/localizations
Npixels = 17; %size of the ROI
bg = 10; %bg photons per pixel
dx = 132; %nm pixel size
dy = 132; %nm
dz=cal.cspline.dz;  %coordinate system of spline PSF is corner based and in units pixels / planes
z0=cal.cspline.z0; % distance and midpoint of stack in spline PSF, needed to translate into nm coordinates
sCMOSvarmap=0; %if scalar : use EMCCD fitter;
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
       xlim([ztruth(1) ztruth(end)])
    subplot(1,2,2)
    plot(ztruth, fractionmisassigned,'-o')
    title('fraction of misassigned localizations')
    xlabel('z (nm)')
    ylabel('fraction')
 
