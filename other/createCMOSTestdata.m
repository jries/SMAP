close all
clearvars

% parameters
    xsize = 32;
    ysize = xsize;
    
    n_frames = 10;
    
    % camera  
        val_offset = 1000; % in counts
        val_noise = 6.5; % in electrons
        val_epc = 0.5; % electrons per count = 1/gain
        val_dark_current = 1; % electrons per frame
        
    % gaussian
        x0 = [xsize/4] .* [2];
        y0 = [xsize/4] .* [2];
        photons = 1000;
        sigma = 1.2; % in pixels
        
rng(1); % reset the random number generator for reproduceable results
    

% maps
    map_pattern_noise = zeros(xsize, ysize); % start with zeros
    map_pattern_noise(eye(size(map_pattern_noise)) == 1) = 1; % create diagonal being white
    map_pattern_noise(end/2:end,:) = 1; % % create 1/2th of area being white
    map_pattern_noise(end/2:2:end,1:2:end/2) = 0; % % create 1/4th of area being checkerboard 
% map_pattern_noise = normrnd(1, 0.1*ones(xsize, ysize));
% map_pattern_noise(17,17) = 23/6.5; % pixel with high readnoise
    
%     map_pattern_offset = zeros(xsize, ysize); % start with zeros     
%     map_pattern_offset(transpose(eye(size(map_pattern_offset))) == 1) = 1; % create diagonal being white
%     map_pattern_offset = repmat(map_pattern_offset(1:end/2, 1:end/2), 2, 2); %
%     map_pattern_offset((end/2+1):2:end,2:2:end/2) = 1; % % create 1/4th of area being checkerboard 
%     map_pattern_offset = fliplr(map_pattern_offset) + 1;
map_pattern_offset = normrnd(1, 0.05*ones(xsize,ysize));
    
    map_pattern_gain = zeros(xsize, ysize); % start with zeros
    map_pattern_gain(end/4:3*end/4, end/4:3*end/4) = 1; % make white center
    map_pattern_gain(transpose(eye(size(map_pattern_gain))) == 1) = 0; % create diagonal being black
    map_pattern_gain(1:2:end/2,1:2:end/2) = 1; % % create 1/4th of area being checkerboard 
map_pattern_gain = normrnd(0, 0.05*ones(xsize, ysize));
    
map_pattern_dark_current =normrnd(1, 0.01*ones(xsize, ysize));
map_pattern_dark_current(17,15) = 500; % pixel with high dark current
    
    
    map_offset = map_pattern_offset * val_offset;
    map_noise = map_pattern_noise * val_noise;
    map_gain = (1 + map_pattern_gain) * val_epc;
    map_dark_current = map_pattern_dark_current * val_dark_current;
    
    % place emitters
        map_gauss = zeros(xsize, ysize); % preallocation
        for xc = x0
            for yc = y0
                for x = 1:xsize
                    for y = 1:ysize
                        map_gauss(x,y) = map_gauss(x,y) + 1/(2*pi*sigma^2) * exp(-(x-xc)^2/2/sigma^2-(y-yc)^2/2/sigma^2);
                    end % for y
                end % for x
            end % for y0
        end % for x0
        
% simulate signal
    for frame = 1:2000
         map_signal(:,:,frame) = map_offset + 1./map_gain .* ( normrnd(0, map_noise) + poissrnd(map_dark_current) + photons*map_gauss );
         map_signal_pure(:,:,frame) = photons*map_gauss;
         map_signal_imphot(:,:,frame) = (map_signal(:,:,frame)-val_offset)*val_epc;
    end % for
    map_signal(map_signal < 0) = 0;
    
    frame = 1;
 
% test plots
    figure
        subplot(3,2,1)
        imagesc(map_gain)
        axis image
        colormap gray
        colorbar
        title('Gain map')
        
        subplot(3,2,2)
        imagesc(map_offset + map_dark_current)
        axis image
        colormap gray
        colorbar
        title('Offset map')
        
        subplot(3,2,3)
        imagesc(map_noise.^2+1)
        axis image
        colormap gray
        colorbar
        title('Variance map')
        
        subplot(3,2,4)
        imagesc(map_dark_current)
        axis image
        colormap jet
        colorbar
        title('Dark current map')
        
        
        subplot(3,2,5)
        imagesc(map_signal(:,:,frame))
        axis image
        colormap jet
        colorbar
        title('Simulated signal')
        
        subplot(3,2,6)
        imagesc(map_signal_pure(:,:,frame))
        axis image
        colormap jet
        colorbar
        title('Photons')
    
    set(gcf, 'Position', [1450 1110 1111 1135])
        return
% write files
    current_dir = cd;
    cd('/Users/diekmann/Documents')
        
    
    tiffstack_filename = '00_CMOS_maps_2HotPX.tif';
    
    % first frame, overwrite the existing file
        imwrite(uint16(map_signal(:, :, 1)), ...
            tiffstack_filename, ...
            'WriteMode', 'overwrite', ...
            'Compression', 'none');
        pause(0.3)
    for n=2:size(map_signal,3)
        imwrite(uint16(map_signal(:, :, n)), ...
                tiffstack_filename, ...
                'WriteMode', 'append', ...
                'Compression', 'none');
            n
            pause(0.001)
    end
    
    gainmap = map_gain;
    offsetmap = map_offset + map_dark_current;
    varmap = map_dark_current + map_noise.^2 +1;
    save('00_CMOS_maps_2HotPX.mat', 'gainmap', 'offsetmap', 'varmap')
    
    
    cd(current_dir);