%% Generating super-resolution images
%%
%% Input:  Nx               : length of image in x axis (px)  
%%         Ny               : length of image in x axis (px)
%%         est_photon       : estimated photons 
%%         est_pos          : estimated positions 
%%         up_ren           : upsampling factor for SR images
%%         width            : FWHM of Gaussian functions for rendering
%%
%% Output: sr_img           : super-resolution image rendered by normalized Gaussian functions  
%%         sr_img_scale     : super-resolution image rendered by weighted Gaussian functions by estimated photons 
%%         sr_img_bin       : super-resolution image as 2D histogram

function [sr_img,sr_img_scale,sr_img_bin] =  super_render(Nx,Ny,est_photon,est_pos,up_ren,width)

%% initailize
sr_img= zeros(Nx*up_ren,Ny*up_ren);
sr_img_scale = zeros(Nx*up_ren,Ny*up_ren);
kernel_size = 3;
Gsigma_sr = width/up_ren/2.3458;     
ss= 1/(2*up_ren);
ind_sr = linspace(ss,kernel_size-ss,kernel_size*up_ren);             
if mod(length(ind_sr),2) == 0
    ind_sr = ind_sr(1:end-1);
end
s = floor(length(ind_sr)/2+1); 
ind_c = mean(ind_sr);

% mesggrid for super-resolution image
[x_grid_sr, y_grid_sr] = meshgrid(ind_sr,ind_sr);                  
% discard locations near boundary
ind = find( (est_pos(:,1)>ceil(kernel_size/2)).*(est_pos(:,2)>ceil(kernel_size/2))...
            .*(est_pos(:,1)< (Ny-ceil(kernel_size/2))).* (est_pos(:,2)< (Nx-ceil(kernel_size/2))));
est_pos = est_pos(ind,:);
est_photon = est_photon(ind);
%%  generate low_res image and high_res images

for ii = 1:length(est_photon)
    pos_bin = round(est_pos(ii,:)*up_ren);
    pos_delta = (est_pos(ii,:)*up_ren-pos_bin)/up_ren;    

    high_kernel = Gauss_kernel(x_grid_sr,y_grid_sr,ind_c+pos_delta(1),ind_c+pos_delta(2),Gsigma_sr);
    high_kernel = high_kernel/sum(high_kernel(:));

    sr_img(pos_bin(2)-s+(1:length(ind_sr)),pos_bin(1)-s+(1:length(ind_sr))) = ...
    sr_img(pos_bin(2)-s+(1:length(ind_sr)),pos_bin(1)-s+(1:length(ind_sr))) + high_kernel;

    sr_img_scale(pos_bin(2)-s+(1:length(ind_sr)),pos_bin(1)-s+(1:length(ind_sr))) = ...
    sr_img_scale(pos_bin(2)-s+(1:length(ind_sr)),pos_bin(1)-s+(1:length(ind_sr))) + est_photon(ii)*high_kernel;
    
end

xLims = linspace(0,Nx*up_ren,Nx*up_ren);    
yLims = linspace(0,Ny*up_ren,Ny*up_ren);    
sr_img_bin = hist3(est_pos*up_ren-0.5,{yLims, xLims});
sr_img_bin = rot90(flipud(sr_img_bin),3);
end