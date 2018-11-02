%% ****************************************************************************************************
%%
%% Fast localization algorithm based on a continuous-space formulation (FALCON) for high density super-resolution microscopy
%% This code was demonstrated in MATLAB 2013a version with GPU
%%
%% Input:  filename          : file name of raw camera images, ex) xxx.tif
%%         numFrame          : number of frames to be reconstructed
%%         dummyFrame        : number of frames to be discarded first
%%         ADU               : photons per camera unit (important)
%%         baseline          : baseline of camera, ex) 100
%%         pixel_size        : raw camera pixel size in nm
%%         EM                : ON = 1 OFF = 0
%%         Gsigma            : widths of PSF (single Gaussian model)
%%         Speed             : option for reconstruction speed
%%         debug             : debug mode
%%
%% Output: Results           : || frame number || x positions || y positions || photon counts || PSF_width_ratio
%%         avgImg            : Temporally averaged image
%% ****************************************************************************************************
function [Results,avgImg]= FALCON_GPU_rel3_WF(imagephot,Gsigma,speed,debug)
%% initialization
dummyFrame=0;
ADU=1;
baseline=0;
EM=0;

beep off;
if isempty(gcp('nocreate')) 
    parpool; 
end

[y_dim,x_dim]= size(imagephot);
%% Common parameters settings
% parameters are opimized for the condition: PSF_width/pixelsize >> 2
boundary = 4;
frame_subset_length = 20;
wavelet_level = 5;                                                         % for background estimation
thresh_delta = 1.1;                                                        % maximum displacement of PSF
thresh_level = 0.05;                                                       % Threshold level for initial localization

% sparsity level (if too many false positive, set para = 2.5 or 3) 
para = 1.5;    

% 
% if EM > 0
%     ADU = ADU/1.4;
% end
if strcmp(speed,'fast')|(speed == 1)
    
    up_decon = 2;
    up_refine = 3;
    iter_step1 = 300;
    iter_step2 = 80;
    iter_step3 = 25;
    bg_iter1 = 50;
    bg_iter2 = 40;
    bg_s = 150;
    
elseif strcmp(speed,'normal')|(speed == 2)
    
    up_decon = 3;
    up_refine = 3;
    iter_step1 = 500;
    iter_step2 = 140;
    iter_step3 = 20;
    bg_iter1 = 100;
    bg_iter2 = 70;
    bg_s = 300;
    
else
    
    up_decon = 3;
    up_refine = 5;
    iter_step1 = 600;
    iter_step2 = 140;
    iter_step3 = 30;
    bg_iter1 = 100;
    bg_iter2 = 70;
    bg_s = 300;
    
end

%% reconstruction grid size

est_y_dim_decon = y_dim*up_decon;
est_x_dim_decon = x_dim*up_decon;
est_y_dim_decon_up = (y_dim)*up_decon;
est_x_dim_decon_up = (x_dim)*up_decon;
est_y_dim_refine = y_dim*up_refine;
est_x_dim_refine = x_dim*up_refine;
est_y_dim_refine_up = (y_dim)*up_refine;
est_x_dim_refine_up = (x_dim)*up_refine;


%% PSF kernels generation
sigma_delta = 0.15;
[fPSF_decon,fPSF_refine,fPSF_dev_x,fPSF_dev_y,fPSF_dev_z] = ...
    Gen_kernelsVer2(x_dim,y_dim,up_decon,up_refine,Gsigma,sigma_delta);

% Fourier coefficients
InvPSF_base = gpuArray(repmat(abs(fPSF_decon).^2,[1,1,frame_subset_length]))/(up_decon)^2;
fPSF_decon_base = gpuArray(repmat(fPSF_decon,[1,1,frame_subset_length]));
fPSF_refine_base = gpuArray(repmat(fPSF_refine,[1,1,frame_subset_length]));
fPSF_dev_x_base = gpuArray(repmat(fPSF_dev_x,[1,1,frame_subset_length]));
fPSF_dev_y_base = gpuArray(repmat(fPSF_dev_y,[1,1,frame_subset_length]));
fPSF_dev_z_base = gpuArray(repmat(fPSF_dev_z,[1,1,frame_subset_length]));

%% downsampling index
% downsampling index for deconvolution
ind_sam = round(up_decon/2);
ind_down_y1 = ind_sam:up_decon:y_dim*up_decon;
ind_down_x1 = ind_sam:up_decon:x_dim*up_decon;

% downsampling index for refinement
ind_sam = round(up_refine/2);
ind_down_y2 = ind_sam:up_refine:y_dim*up_refine;
ind_down_x2 = ind_sam:up_refine:x_dim*up_refine;

% macro function
D1 = @(x) x(ind_down_y1,ind_down_x1,:);
D2 = @(x) x(ind_down_y2,ind_down_x2,:);

%% output
numFrame=1;
Results = zeros(1000*numFrame,5,'single');
index_Results = 1;
avgImg = zeros(y_dim,x_dim,'single');

%% localization start
% fprintf('FALCON GPU start...\n');
for frame_count = 1: ceil(numFrame/frame_subset_length)
    
    sframe = (frame_count-1)*frame_subset_length+1+dummyFrame;
    eframe = min(frame_count*frame_subset_length,numFrame)+dummyFrame;
    subset_length = eframe-sframe+1;
    CCD_imgs_sub = zeros(y_dim,x_dim,subset_length,'single');
    
%     fprintf('frames no. %d ~ %d are reconstructing...\n',sframe,eframe);
    
    % load camera images
    for zz = sframe:eframe
         CCD_imgs_sub(:,:,zz-sframe+1) = imagephot(:,:,zz);
    end
    % averaged image
    avgImg = avgImg +sum(CCD_imgs_sub,3)*ADU/numFrame;
    % convert camera unit in photon counts
    CCD_imgs_sub = (max(CCD_imgs_sub-baseline,0))*ADU;
    % initial background estimation
    for tt= 1:subset_length
        medVal = median(CCD_imgs_sub(:))*2;
        CCD_imgs_sub_temp = CCD_imgs_sub;
        ind = CCD_imgs_sub_temp > medVal;
        CCD_imgs_sub_temp(ind) = medVal;
    end
    est_backgrounds = mywaveletfilter(CCD_imgs_sub_temp,4); %XXX adjust wavelet level?
    
    
    % initialize variables
    est_c = gpuArray(zeros(est_y_dim_decon_up,est_x_dim_decon_up,subset_length,'single'));          % super resolution grid.
    d = gpuArray(zeros(est_y_dim_decon_up,est_x_dim_decon_up,subset_length,'single'));              % auxilliary variable
    v = gpuArray(zeros(est_y_dim_decon_up,est_x_dim_decon_up,subset_length,'single'));              % auxilliary variable
    gradient_c = gpuArray(zeros(est_y_dim_refine_up,est_x_dim_refine_up,subset_length,'single'));
    gradient_dx = gpuArray(zeros(est_y_dim_refine_up,est_x_dim_refine_up,subset_length,'single'));
    gradient_dy = gpuArray(zeros(est_y_dim_refine_up,est_x_dim_refine_up,subset_length,'single'));
    gradient_dz = gpuArray(zeros(est_y_dim_refine_up,est_x_dim_refine_up,subset_length,'single'));
    grad_c_temp = gpuArray(zeros(est_y_dim_refine_up,est_x_dim_refine_up,subset_length,'single'));
    grad_d_temp = gpuArray(zeros(est_y_dim_refine_up,est_x_dim_refine_up,subset_length,'single'));
    
    % parameter settings for ADMM
    bg_min = max(CCD_imgs_sub(:))*0.02;
    w = para*sqrt(max(imresize(est_backgrounds,up_decon),bg_min));         % weights for sparsity priors
    mu = 0.05;                                                             % auxiliary parameter for ADMM
    thresh_c = gpuArray(w/mu);
    
    % PSF kernels in Fourier domain.
    InvPSF = 1./(InvPSF_base(:,:,1:subset_length) + mu*1);
    fPSF_decon = fPSF_decon_base(:,:,1:subset_length);
    fPSF_refine = fPSF_refine_base(:,:,1:subset_length);
    fPSF_dev_x = fPSF_dev_x_base(:,:,1:subset_length);
    fPSF_dev_y = fPSF_dev_y_base(:,:,1:subset_length);
    fPSF_dev_z = fPSF_dev_z_base(:,:,1:subset_length);
    fPSF_dev_x_t = -fPSF_dev_x;
    fPSF_dev_y_t = -fPSF_dev_y;
    fPSF_dev_z_t = fPSF_dev_z;
    
    
    % define macro function
    InvH = @(x) real(ifft2( InvPSF.*fft2(x)));
    A = @(x,delta_x,delta_y,delta_z) real(D2(ifft2( fft2(x).*fPSF_refine -fft2(delta_x.*x).*fPSF_dev_x ...
        - fft2(delta_y.*x).*fPSF_dev_y - fft2(delta_z.*x).*fPSF_dev_z)));
    At = @(fx,delta_x,delta_y,delta_z) real(ifft2(fx.*fPSF_refine) - delta_x.*ifft2(fx.*fPSF_dev_x_t) ...
        - delta_y.*ifft2(fx.*fPSF_dev_y_t)  - delta_z.*ifft2(fx.*fPSF_dev_z_t));
    
    Ad = @(est_c,delta_x,delta_y,delta_z) real(D2(ifft2(fft2(delta_x.*est_c).*fPSF_dev_x ...
        + fft2(delta_y.*est_c).*fPSF_dev_y  + fft2(delta_z.*est_c).*fPSF_dev_z - fft2(est_c).*fPSF_refine )));
    Adx = @(x) real(D2(ifft2(fft2(x).*fPSF_dev_x)));
    Ady = @(x) real(D2(ifft2(fft2(x).*fPSF_dev_y)));
    Adz = @(x) real(D2(ifft2(fft2(x).*fPSF_dev_z)));
    Adxt = @(x) real(ifft2(x.*fPSF_dev_x_t));
    Adyt = @(x) real(ifft2(x.*fPSF_dev_y_t));
    Adzt = @(x) real(ifft2(x.*fPSF_dev_z_t));
    
    %% Step1 Deconvolution with sparse priors
    
    % pre-calculation for the gradients
    CCD_imgs_sub = gpuArray(CCD_imgs_sub);
    est_backgrounds = gpuArray(est_backgrounds);
    CCD_imgs_sub_bg = gpuArray(zeros(est_y_dim_decon_up,est_x_dim_decon_up,subset_length,'single'));
    CCD_imgs_sub_bg(ind_down_y1,ind_down_x1,:) = CCD_imgs_sub-est_backgrounds;
    temp_PSF_t = real(ifft2(fft2(CCD_imgs_sub_bg).*fPSF_decon));
    
    for count = 1:iter_step1
        r =   v + d;
        est_c = InvH(temp_PSF_t+mu*r);
        v_temp = est_c - d;
        v (1:est_y_dim_decon,1:est_x_dim_decon,:) = max(v_temp(1:est_y_dim_decon,1:est_x_dim_decon,:)-thresh_c,0);
        d = d + v- est_c;
        
        % background update
        if (mod(count,bg_iter1) == 0)&&(count>=bg_s)
            res =  gather(CCD_imgs_sub-D1(real(ifft2(fft2(est_c).*fPSF_decon))));
%             est_backgrounds_new = padarray(background_estimation(res((boundary+1):end-boundary,(boundary+1):end-boundary,:),1,wavelet_level,'db6',2),[boundary,boundary],'replicate','both');
             est_backgrounds_new = padarray(mywaveletfilter(res((boundary+1):end-boundary,(boundary+1):end-boundary,:),4),[boundary,boundary],'replicate','both');
            w = para*sqrt(max(imresize(single(est_backgrounds_new),up_decon),bg_min));          % weights for sparsity priors
            thresh_c = gpuArray(w/mu);
            est_backgrounds = gpuArray(est_backgrounds_new);
            CCD_imgs_sub_bg(ind_down_y1,ind_down_x1,:) = CCD_imgs_sub-est_backgrounds;
            temp_PSF_t = real(ifft2(fft2(CCD_imgs_sub_bg).*fPSF_decon));
        end
        
        % debug mode : display current estimated variables
        if debug > 0
            figure(100);
            subplot(1,3,1,'align')
            imagesc(CCD_imgs_sub(:,:,1)); colormap(hot); axis image; axis off; axis tight; colormap(hot); title('CCD img'); colorbar;
            subplot(1,3,2,'align')
            imagesc(est_c((boundary*up_decon+1):end-boundary*up_decon,(boundary*up_decon+1):end-boundary*up_decon,1),[0 200]);
            colormap(hot); axis image; axis off; axis tight; title(['SR img, iter = ',num2str(count)]); colorbar;
            subplot(1,3,3,'align')
            imagesc(est_backgrounds((boundary+1):end-boundary,(boundary+1):end-boundary,1));
            colormap(hot); axis image; axis off; axis tight; title('current background'); colorbar;
        end
    end
    
    %% Step2 Deconvolution with spatial support
    thresh = mean(max(max(est_c,[],1),[],2))*thresh_level;
    est_c(est_c<thresh) = 0;
    ind_set = est_c>0;
    v(~ind_set) = 0;
    if sum(ind_set(:))>0
        % update est_c on the support.
        for count = 1:iter_step2
            r =   v + d;
            est_c = InvH(temp_PSF_t+mu*r);
            v(ind_set) = max(est_c(ind_set) - d(ind_set),0);
            d = d + v - est_c;
            
            % background update
            if mod(count,bg_iter2) == 0
                res =  gather(CCD_imgs_sub-D1(real(ifft2(fft2(est_c).*fPSF_decon))));
%                 est_backgrounds_new = padarray(background_estimation(res((boundary+1):end-boundary,(boundary+1):end-boundary,:),0,wavelet_level,'db6',1),[boundary,boundary],'replicate','both');
                est_backgrounds_new = padarray(mywaveletfilter(res((boundary+1):end-boundary,(boundary+1):end-boundary,:),4),[boundary,boundary],'replicate','both');
                est_backgrounds = gpuArray(est_backgrounds_new);
                CCD_imgs_sub_bg(ind_down_y1,ind_down_x1,:) = CCD_imgs_sub-est_backgrounds;
                temp_PSF_t = real(ifft2(fft2(CCD_imgs_sub_bg).*fPSF_decon));
            end
            
            % debug mode : display current estimated variables
            if debug > 0
                figure(101);
                subplot(1,3,1,'align')
                imagesc(CCD_imgs_sub(:,:,1)); colormap(hot); axis image; axis off; axis tight; colormap(hot); title('CCD img'); colorbar;
                subplot(1,3,2,'align')
                imagesc(est_c((boundary*up_decon+1):end-boundary*up_decon,(boundary*up_decon+1):end-boundary*up_decon,1),[0 200]);
                colormap(hot); axis image; axis off; axis tight; title(['SR img, iter = ',num2str(count)]); colorbar;
                subplot(1,3,3,'align')
                imagesc(est_backgrounds((boundary+1):end-boundary,(boundary+1):end-boundary,1));
                colormap(hot); axis image; axis off; axis tight; title('current background'); colorbar;
            end
        end
        
        
        %% Step3 Continuous refinement
        
        % get initial localizations
        thresh = mean(max(max(est_c,[],1),[],2))*thresh_level;
        est_c(est_c<0) = 0;
        [est_c,delta_x,delta_y] = Peakdets(gather(est_c),thresh);
        ind_set = find(est_c>0);
        
        % grid rearrangement
        if up_decon ~= up_refine
            [est_pos_y,est_pos_x,nn] = ind2sub(size(est_c),ind_set);
            est_photon = est_c(ind_set);
            offset1 = mod(up_decon/2,1);
            offset2 = mod(up_refine/2,1);
            est_pos_x = (est_pos_x+delta_x(ind_set)-offset1)/up_decon*up_refine+offset2;
            est_pos_y = (est_pos_y+delta_y(ind_set)-offset1)/up_decon*up_refine+offset2;
            ind_set = sub2ind([est_y_dim_refine_up, est_x_dim_refine_up,subset_length],round(est_pos_y),round(est_pos_x),nn);
            est_c = gpuArray(zeros(est_y_dim_refine_up, est_x_dim_refine_up,subset_length,'single'));
            delta_x = gpuArray(zeros(est_y_dim_refine_up, est_x_dim_refine_up,subset_length,'single'));
            delta_y = gpuArray(zeros(est_y_dim_refine_up, est_x_dim_refine_up,subset_length,'single'));
            est_c(ind_set) = (est_photon);
            delta_x(ind_set) = (est_pos_x-round(est_pos_x));
            delta_y(ind_set) = (est_pos_y-round(est_pos_y));
            CCD_imgs_sub_bg = gpuArray(zeros(est_y_dim_refine_up,est_x_dim_refine_up,subset_length,'single'));
            CCD_imgs_sub_bg(ind_down_y2,ind_down_x2,:) = CCD_imgs_sub-est_backgrounds;
        end
        
        % variables to be estimated
        est_c = gpuArray(est_c);
        delta_x = gpuArray(delta_x);
        delta_y = gpuArray(delta_y);
        delta_z = gpuArray(zeros(est_y_dim_refine_up, est_x_dim_refine_up,subset_length,'single'));
        step_boost = 1.6;
        
        % pre-calculation for the gradients
        Id = @(x) x(ind_set);
        temp_PSF_t = (real(ifft2(fft2(CCD_imgs_sub_bg).*fPSF_refine)));
        temp_PSF_dev_x_t = (real(ifft2(fft2(CCD_imgs_sub_bg).*fPSF_dev_x_t)));
        temp_PSF_dev_y_t = (real(ifft2(fft2(CCD_imgs_sub_bg).*fPSF_dev_y_t)));
        temp_PSF_dev_z_t = (real(ifft2(fft2(CCD_imgs_sub_bg).*fPSF_dev_z_t)));
        
        
        
        for count = 1:iter_step3
            
            % debug mode : display current estimated variables
            if debug > 0
                figure(102);
                subplot(2,2,1,'align')
                imagesc(est_c((boundary*up_decon+1):end-boundary*up_decon,(boundary*up_decon+1):end-boundary*up_decon,1));
                colormap(hot); axis image; axis off; axis tight; title(['SR img, iter = ',num2str(count)]); colorbar;
                subplot(2,2,2,'align')
                imagesc(delta_x((boundary*up_decon+1):end-boundary*up_decon,(boundary*up_decon+1):end-boundary*up_decon,1));
                colormap(hot); axis image; axis off; axis tight; title(['delta x, iter = ',num2str(count)]); colorbar;
                subplot(2,2,3,'align')
                imagesc(delta_y((boundary*up_decon+1):end-boundary*up_decon,(boundary*up_decon+1):end-boundary*up_decon,1));
                colormap(hot); axis image; axis off; axis tight; title(['delta y, iter = ',num2str(count)]); colorbar;
                subplot(2,2,4,'align')
                imagesc(delta_z((boundary*up_decon+1):end-boundary*up_decon,(boundary*up_decon+1):end-boundary*up_decon,1));
                colormap(hot); axis image; axis off; axis tight; title(['delta witdh, iter = ',num2str(count)]); colorbar;
            end
            
            % % refine locations in x,y and width of PSF
            grad_d_temp(ind_down_y2,ind_down_x2,:) = Ad(est_c, delta_x,delta_y,delta_z);
            fgrad_d_temp = fft2(grad_d_temp);
            gradient_dx(ind_set) = 2*est_c(ind_set).*(Id(Adxt(fgrad_d_temp))+Id(temp_PSF_dev_x_t));
            gradient_dy(ind_set) = 2*est_c(ind_set).*(Id(Adyt(fgrad_d_temp))+Id(temp_PSF_dev_y_t));
            
            % determine step size
            temp_res =  (D2(CCD_imgs_sub_bg) + D2(grad_d_temp));
            temp_dx = Adx(est_c.*gradient_dx);
            temp_dy = Ady(est_c.*gradient_dy);
            for tt = 1:subset_length
                temp_res_sub =temp_res((boundary+1):end-boundary,(boundary+1):end-boundary,tt);
                temp_dx_sub = temp_dx((boundary+1):end-boundary,(boundary+1):end-boundary,tt);
                temp_dy_sub = temp_dy((boundary+1):end-boundary,(boundary+1):end-boundary,tt);
                gradient_dx(:,:,tt) = gradient_dx(:,:,tt).*max((temp_res_sub(:)'*temp_dx_sub(:))/(temp_dx_sub(:)'*temp_dx_sub(:)),0);
                gradient_dy(:,:,tt) = gradient_dy(:,:,tt).*max((temp_res_sub(:)'*temp_dy_sub(:))/(temp_dy_sub(:)'*temp_dy_sub(:)),0);
            end
            
            % update
            delta_x(ind_set) = delta_x(ind_set) - step_boost*gradient_dx(ind_set);
            delta_y(ind_set) = delta_y(ind_set) - step_boost*gradient_dy(ind_set);
            delta_x(ind_set) = sign(delta_x(ind_set)).*min(abs(delta_x(ind_set)),thresh_delta);
            delta_y(ind_set) = sign(delta_y(ind_set)).*min(abs(delta_y(ind_set)),thresh_delta);
            
            % rearrangement
            [est_c,delta_x,delta_y,delta_z] = handover(est_c,delta_x,delta_y,delta_z);
            ind_set = find(est_c>0);
            Id = @(x) x(ind_set);
            
            
            % % refine photons
            grad_c_temp(ind_down_y2,ind_down_x2,:) = A(est_c,delta_x,delta_y,delta_z);
            gradient_c(ind_set) = Id(At(fft2(grad_c_temp),delta_x,delta_y,delta_z)- temp_PSF_t)  + delta_x(ind_set).*Id(temp_PSF_dev_x_t)...
                + delta_y(ind_set).*Id(temp_PSF_dev_y_t) + delta_z(ind_set).*Id(temp_PSF_dev_z_t);
            % determine step size
            temp_res =  (D2(CCD_imgs_sub_bg) - D2(grad_c_temp));
            temp_c = (A(gradient_c,delta_x,delta_y,delta_z));
            for tt = 1:subset_length
                temp_res_sub = temp_res((boundary+1):end-boundary,(boundary+1):end-boundary,tt);
                temp_c_sub = temp_c((boundary+1):end-boundary,(boundary+1):end-boundary,tt);
                gradient_c(:,:,tt) = gradient_c(:,:,tt).*max(-(temp_res_sub(:)'*temp_c_sub(:))/(temp_c_sub(:)'*temp_c_sub(:)),0);
            end
            % update
            est_c(ind_set) = max(est_c(ind_set) - gradient_c(ind_set),1);
            
            
            
            % % refine width of PSF
            grad_d_temp(ind_down_y2,ind_down_x2,:) = Ad(est_c, delta_x,delta_y,delta_z);
            fgrad_d_temp = fft2(grad_d_temp);
            gradient_dz(ind_set) = 2*est_c(ind_set).*(Id(Adzt(fgrad_d_temp))+Id(temp_PSF_dev_z_t));
            
            % determine step size
            temp_res =  (D2(CCD_imgs_sub_bg) + D2(grad_d_temp));
            temp_dz = Adz(est_c.*gradient_dz);
            for tt = 1:subset_length
                temp_res_sub =temp_res((boundary+1):end-boundary,(boundary+1):end-boundary,tt);
                temp_dz_sub = temp_dz((boundary+1):end-boundary,(boundary+1):end-boundary,tt);
                gradient_dz(:,:,tt) = gradient_dz(:,:,tt).*max((temp_res_sub(:)'*temp_dz_sub(:))/(temp_dz_sub(:)'*temp_dz_sub(:)),0);
            end
            
            % update
            delta_z(ind_set) = delta_z(ind_set) - gradient_dz(ind_set);
            delta_z(ind_set) = sign(delta_z(ind_set)).*(min(abs(delta_z(ind_set)),2*sigma_delta));
            
            
        end
        
    end
    
    
    %% localization results save
    thresh = mean(max(max(est_c,[],1),[],2))*thresh_level;
    est_c(est_c<thresh) = 0;
    mask = zeros(est_y_dim_refine,est_x_dim_refine,subset_length);
    mask(boundary*up_refine+1:end-boundary*up_refine,boundary*up_refine+1:end-boundary*up_refine,:) = 1;
    est_c = gather(est_c(1:est_y_dim_refine,1:est_x_dim_refine,:).*mask);
    delta_x = gather(delta_x(1:est_y_dim_refine,1:est_x_dim_refine,:).*mask);
    delta_y = gather(delta_y(1:est_y_dim_refine,1:est_x_dim_refine,:).*mask);
    delta_z = gather(delta_z(1:est_y_dim_refine,1:est_x_dim_refine,:).*mask);
    
    ind_set = find((est_c>1));
    est_photon = est_c(ind_set);
    if EM > 0
        est_photon = est_photon*1.4;
    end
    est_sigma_delta = 1+delta_z(ind_set);
    [est_pos_y, est_pos_x,numFrame_sub] = ind2sub(size(est_c),ind_set);
    est_pos = [est_pos_x+delta_x(ind_set) est_pos_y+delta_y(ind_set)];
    offset = mod(up_refine/2,1);
    % convert coordinates from sub-pixel to CCD pixel
    est_pos = (est_pos-offset)/up_refine;
    
    num_cluster = length(ind_set);
    ee = index_Results+num_cluster-1;
    if ee > size(Results,1)
        Results = [Results;zeros(10000,5)];
    end
    % save results
    Results(index_Results:ee,1) =  (frame_count-1)*frame_subset_length +numFrame_sub;
    Results(index_Results:ee,2:3) = single(est_pos);
    Results(index_Results:ee,4) = single(est_photon);
    Results(index_Results:ee,5) = single(est_sigma_delta);
    if mod(frame_count,10) == 0
        save('FALCON_temp.mat','Results','avgImg');
    end
    index_Results = ee+1;
end
ind = Results(:,4) > 1;
Results = double(Results(ind,:));
end

