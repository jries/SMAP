%  Copyright (c)2017 Ries Lab, European Molecular Biology Laboratory,
%  Heidelberg.
%  
%  Additional permission under GNU GPL version 3 section 7
%  
%  If you modify this Program, or any covered work, by
%  linking or combining it with libraries required for interaction
%  with analysis programs such as Igor Pro or Matlab,
%  the licensors of this Program grant you additional permission
%  to convey the resulting work.

%% example script for calibration and correction of depth-dependent aberrations
%Installation instructions
%  Requires Matlab 2016a or newer including the following tool boxes:
%  image analysis, curve fitting, optimization, statistics

    
%% depth-dependent calibration
% fit bead stacks in gel as shown above using the same fitter with the
% same, save fitted z-positions and frames
% z-positions should not be corrected for the refractive index mismatch. If
% they are, divide by the RI_mismatch_factor.

%load fitted coordinates of beads in gel
beads=(readtable('data/beads_in_gel_2.csv'));
beads.z=-beads.z; %revert z: in this example beads close to the coverslip have a larger z: this comes from the calibration of the PSF, which uses objective positions
% paramters for bead stack
p.dz=20;
% determine 3d corrections
zcorr=calibrate3Daberrations(beads,p);
save('data/3Daberrationcorrection.mat','zcorr') %save correction to be used later

%% correct for depth-induced errors
load('data/3Daberrationcorrection.mat') %correction
% fit your data with a 3D method of choice (example saved in
% 'data/cme_2um_site_uncorrected_sml.csv'
locs=readtable('data/cme_2um_site_uncorrected_sml.csv');
locs.z=-locs.z;

RI_mismatch_factor=.75; %refractive index mismatch correction

binsize=10;

[img2d,xedges,yedges]=histcounts2(locs.x,locs.y,'BinWidth',[1 1]*binsize); %reconstruct superresolution image
h=fspecial('gauss',5,.6); %and blur a bit for rendering
figure(104);subplot(2,2,1);imagesc(xedges,yedges,filter2(h,img2d'));
axis equal
title('x-y')
xlabel('x (nm)');ylabel('y (nm)');
%%plot side view reconstruction
ypos=250; %position of slice in nm
slicewidth=80; %width of slice in nm
iny=locs.y>ypos-slicewidth/2 & locs.y<ypos+slicewidth/2; 
[img2dz,xedges,zedges]=histcounts2(locs.x(iny),locs.z(iny)*RI_mismatch_factor,'BinWidth',[1 1]*binsize);
figure(104);subplot(2,2,2);imagesc(xedges,zedges,filter2(h,img2dz'));
axis equal
title('x-z raw')
xlabel('x (nm)');ylabel('z (nm)');
% correct for aberrations and refractive index mismatch
objective_depth=2000; %position of objective above coverslip

locs.z_corr=correct_3Daberrations(zcorr,locs.z,objective_depth)*RI_mismatch_factor; %correct for aberrations
[img2dz,xedges,zedges]=histcounts2(locs.x(iny),locs.z_corr(iny),'BinWidth',[1 1]*binsize); %plot sideview
figure(104);subplot(2,2,4);imagesc(xedges,zedges,filter2(h,img2dz'));
title('x-z corrected')
xlabel('x (nm)');ylabel('z (nm)');
colormap gray
axis equal