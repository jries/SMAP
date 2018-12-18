clearvars
close all

% parameters
    xsize = 32;
    ysize = xsize;
    
    n_frames = 10;
    
    % camera  
        val_offset = 100; % in counts
        val_noise = 10; % in electrons
        val_epc = 0.5; % electrons per count = 1/gain
    
    % gaussian
        x0 = [xsize/4] .* [1 3];
        y0 = [xsize/4] .* [3 1];
        photons = 2000;
        sigma = 3; % in pixels
        
rng(1); % reset the random number generator for reproduceable results
    

% maps
    map_pattern_noise = zeros(xsize, ysize); % start with zeros
    map_pattern_noise(eye(size(map_pattern_noise)) == 1) = 1; % create diagonal being white
    map_pattern_noise(end/2:end,:) = 1; % % create 1/2th of area being white
    map_pattern_noise(end/2:2:end,1:2:end/2) = 0; % % create 1/4th of area being checkerboard 
    
    map_pattern_offset = zeros(xsize, ysize); % start with zeros     
    map_pattern_offset(transpose(eye(size(map_pattern_offset))) == 1) = 1; % create diagonal being white
    map_pattern_offset = repmat(map_pattern_offset(1:end/2, 1:end/2), 2, 2); %
    map_pattern_offset((end/2+1):2:end,2:2:end/2) = 1; % % create 1/4th of area being checkerboard 
    map_pattern_offset = fliplr(map_pattern_offset);
    
    map_pattern_gain = zeros(xsize, ysize); % start with zeros
    map_pattern_gain(end/4:3*end/4, end/4:3*end/4) = 1; % make white center
    map_pattern_gain(transpose(eye(size(map_pattern_gain))) == 1) = 0; % create diagonal being black
    map_pattern_gain(1:2:end/2,1:2:end/2) = 1; % % create 1/4th of area being checkerboard 

    map_offset = map_pattern_offset * val_offset;
    map_noise = map_pattern_noise * val_noise;
    map_gain = (1 + map_pattern_gain) / val_epc;
    
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
    for frame = 1
         map_signal(:,:,frame) = map_offset + map_gain .* ( normrnd(0, map_noise) + poissrnd(photons*map_gauss) );
    end % for
 
% test plots
    figure
        imagesc(map_gain)
        axis image
        colormap gray
        colorbar
        
    figure
        imagesc(map_offset)
        axis image
        colormap gray
        colorbar
        
    figure
        imagesc(map_noise)
        axis image
        colormap gray
        colorbar
        
    figure
        imagesc(map_signal(:,:,frame))
        axis image
        colormap gray
        colorbar
        
        
% write files
    current_dir = cd;
    cd('/Users/diekmann/Documents')
        
    imwrite(uint16(map_signal), 'CMOS_testsignal.tif')
    
    cd(current_dir)