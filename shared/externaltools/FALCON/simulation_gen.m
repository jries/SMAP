%% ****************************************************************************************************
%% Simulated image generation
%%
%% Input:  img_size          : size of CCD image
%%         num_frame         : number of images
%%         boundary_offset   : boundary doesn't have any molecule
%%         num_mol           : number of activated molecule in each frame
%%         mol_mu            : mean value of emitted photons
%%         mol_sigma         : standard deviation of emitted photons
%%         Gsigma1 & 2       : widths of Gaussian fuctions for PSF
%%         Gsimga_ratio      : ratio for two Gaussian functions for PSF: PSF = Gsimga_ratio*F_sigma1 + (1-Gsimga_ratio)*F_sigma2
%%         background        : background fluorescent signal
%%         baseline          : uniform offset for final image.
%%         EM                : EM-gain option 
%%
%% Output: CCD_imgs          : blurred low resolution image + noise [short exposure multiple shots]
%%         true_pos          : 2D coordinate of molecules' position (in px)
%%         true_photon       : Number of emitted photons
%%
%% ****************************************************************************************************

function [CCD_imgs,true_pos,true_photon] =  simulation_gen(img_size,num_frame,boundary_offset,num_mol,mol_mu,mol_sigma,Gsigma1,Gsigma2,Gsigma_ratio,background,readout_rms,baseline,EM)

% CCD images
CCD_imgs = zeros(img_size,img_size,num_frame);                             
% photon emission stat followed by log-normal distribution.

m = log((mol_mu^2)/sqrt(mol_sigma^2+mol_mu^2));                            % mean for log-normal function
s = sqrt(log(mol_sigma^2/(mol_mu^2)+1));                                   % sigma for log-normal function
true_photon = lognrnd(m,s,num_mol,num_frame);                              
true_pos = zeros(num_mol,2,num_frame);

% low-resolution grid
ind_low = linspace(0.5,img_size-0.5,img_size);                             
[x_grid_low, y_grid_low] = meshgrid(ind_low,ind_low);    
% boundary 
true_s = ind_low(boundary_offset)+0.5;
true_e = ind_low(end-boundary_offset)+0.5;


for t = 1:num_frame
    %  random distribution
    true_pos(:,1,t) = true_s+(true_e-true_s).*rand(num_mol,1);               
    true_pos(:,2,t) = true_s+(true_e-true_s).*rand(num_mol,1);
      
%     %  Synthetic ring 
%     radius = 2;
%     theta = 2*pi*rand(num_mol,1);
%     true_pos(:,1,t) = img_size/2 + radius*cos(theta);
%     true_pos(:,2,t) = img_size/2 + radius*sin(theta);

    %  generate CCD image and SR images
    for ii = 1:num_mol
        low_kernel = Gsigma_ratio*Gauss_kernel(x_grid_low,y_grid_low,true_pos(ii,1,t),true_pos(ii,2,t),Gsigma1)...
                     +(1-Gsigma_ratio)*Gauss_kernel(x_grid_low,y_grid_low,true_pos(ii,1,t),true_pos(ii,2,t),Gsigma2);
        low_kernel = low_kernel/sum(low_kernel(:));
        CCD_imgs(:,:,t) = CCD_imgs(:,:,t) + true_photon(ii,t)*low_kernel;
    end
end

CCD_imgs = CCD_imgs + background;                                          % add background fluorescent light
if EM>0
    CCD_imgs = CCD_imgs + 1.4*(poissrnd(CCD_imgs)-CCD_imgs);               % shot-noise with excess noise factor of EM gain
else
    CCD_imgs = poissrnd(CCD_imgs);                                         % shot-noise with excess noise without EM gain
end
CCD_imgs = CCD_imgs + readout_rms*randn(size(CCD_imgs));                   % readout noise followed by Gaussian stat
if EM>0
     CCD_imgs = CCD_imgs*EM;
end 
CCD_imgs = CCD_imgs + baseline;                                            % add offset baseline
CCD_imgs = round(CCD_imgs);                                                % discretization
CCD_imgs(CCD_imgs<0) = 0;
end