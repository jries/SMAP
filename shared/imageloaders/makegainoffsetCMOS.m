function [gainmap,offsetmap,varmap,roi]=makegainoffsetCMOS(camfname,exposuretime_data)
    % if the camera has been characterized exposure time dependent,
    % find behavior for exposure time used here
    if length(camfname)>4 && strcmp(camfname(end-3:end),'.mat')
        camfname=camfname(1:end-4);
    end
    if length(camfname)>5 && strcmp(camfname(end-4:end),'.calb')
        camfname=camfname(1:end-5);
    end
    camfilemat=['settings' filesep 'cameras' filesep camfname '.mat'];
    camfilejson=['settings' filesep 'cameras' filesep camfname '.calb'];
    if exist(camfilemat,'file')
        l=load(camfilemat);
        roi=[];
        if isfield(l, 'read_noise_variance') ...
                && isfield(l, 'thermal_noise_variance_per_s') ...
                && isfield(l, 'pixel_baseline') ...
                && isfield(l, 'thermal_counts_per_s')
            [offsetmap, varmap] = makeExpDependMap(l, exposuretime_data);
            gainmap=1./l.gainmap;
        elseif isfield(l, 'offsetmap')
            offsetmap=l.offsetmap;
            varmap=l.varmap;
            gainmap=1./l.gainmap;
        elseif isfield(l, 'mean')
            offsetmap=l.mean;
            varmap=l.variance;  
            gainmap=1./l.metadata.pix2phot*ones(size(offsetmap));
            roi=l.metadata.roi;
        else 
            disp(['no camera calibration found in file ' camfname])
        end
        %assume gainmap is independent of exposure time
        %NB: gainmap has to be inverted, since we need the unit
        %    electrons/counts
%             gainmap=1/l.gainmap;
        %only use media of gainmap
            gainmap = ones(size(gainmap))*median(gainmap(:));
        %varmap has to be in units of photons^2
            varmap=varmap.*gainmap.^2; 
            
                %???? in units of photons, to be used later directly. 
                %A: Yes, conversion to photons^2 is not considered in GPUmleFIT
    elseif exist(camfilejson,'file')
        txt=fileread(camfilejson);
        sr=jsondecode(txt);
        l.read_noise_variance=reshape(sr.rnSq,sr.width,sr.height)';
        l.thermal_noise_variance_per_s=reshape(sr.tnSqPerSec,sr.width,sr.height)';
        l.pixel_baseline=reshape(sr.baseline,sr.width,sr.height)';
        l.thermal_counts_per_s=reshape(sr.dcPerSec,sr.width,sr.height)';
        [offsetmap, varmap] = makeExpDependMap(l, exposuretime_data);
        %what do we do with gain map? Corrrect? Already averaged?
        gainmap=reshape(sr.gain,sr.width,sr.height)';
            %assume gainmap is independent of exposure time
            %NB: gainmap has to be inverted, since we need the unit
            %    electrons/counts
                gainmap=1./gainmap;
            %only use media of gainmap
                gainmap = ones(size(gainmap))*median(gainmap(:));
        varmap=varmap.*gainmap.^2;
        if isfield(sr,'x') && sr.x>=0 %roi defined 
            roi=[sr.x sr.y sr.width sr.height]; %1 based
        else
            roi=[];
        end
    else
        gainmap=[];offsetmap=[];varmap=[];roi=[];
    end
    
end

function [offsetmap, varmap] = makeExpDependMap(l, exposuretime_data)
        exposure_time_s = exposuretime_data /1e3; % maps have units seconds, SMAP handles timing in ms
        
        varmap = l.read_noise_variance + l.thermal_noise_variance_per_s * exposure_time_s;
        offsetmap = l.pixel_baseline + l.thermal_counts_per_s * exposure_time_s;
end % function