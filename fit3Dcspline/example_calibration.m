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

%% example generating a PSF model from a bead stack
%Installation instructions
%  Requires Matlab 2016a or newer including the following tool boxes:
%  image analysis, curve fitting, optimization, statistics
%  GPU support for CUDA 8

%  The calibration can either be done using the following script or using
%  GUI by running calibrate3D_GUI.m

%% add path to helper functions
addpath('shared')

%% make bead calibration for 3D astigmatic PSF
datapath=[pwd filesep 'example_data' filesep 'beadstacks_3D_astig'];
p.filelist={[datapath filesep 'stack3D_1.tif'],[datapath filesep 'stack3D_2.tif'],[datapath filesep 'stack3D_3.tif']};  % list of image stacks
p.outputfile = [pwd filesep 'example_data' filesep 'bead_astig_3dcal.mat' ]; %output file name
p.filtersize =2; %size of the filter (pixels) for bead segmentation. Use larger value for bi-lobed PSF.
p.mindistance=25; %minimum distance (pixels) between beads
p.dz =10; %distance between frames in nm
p.modality ='astigmatism'; 
p.zcorr ='cross-correlation'; %modality of z-alignment
p.zcorrframes=50; % number of frames used for 3D correlation
p.ROIxy = 31; %size of the PSF model in x,y (pixels)

p.smoothz =1; %smoothing parameter in Z
p.gaussrange =[-700 700]; %z range (nm) in which to calibrate the parameters for Gaussian z-fit.
p.gaussroi =19; %size of the ROI in pixels for the Gaussian calibration

calibrate3D(p) %call calibraion function

%% make bead calibration for unmodified PSF
%as before, but using the bead stack with an unmodified PSF
datapath=[pwd filesep 'example_data' filesep 'beadstacks_2D'];
p.filelist={[datapath filesep 'stack2D_1.tif'],[datapath filesep 'stack2D_2.tif'],[datapath filesep 'stack2D_3.tif']};  % list of image stacks
p.outputfile = [pwd filesep 'example_data' filesep 'bead_2DPSF_3dcal.mat' ]; %output file name
p.modality ='2D PSF'; 
calibrate3D(p) %call calibraion function