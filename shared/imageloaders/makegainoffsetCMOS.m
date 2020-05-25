function [gainmap,offsetmap,varmap]=makegainoffsetCMOS(camfname,exposuretime_data)
    % if the camera has been characterized exposure time dependent,
    % find behavior for exposure time used here
    if length(camfname)>4 && strcmp(camfname(end-3:end),'.mat')
        camfname=camfname(1:end-4);
    end
    camfile=['settings' filesep 'cameras' filesep camfname '.mat'];
    if exist(camfile,'file')
        l=load(camfile);
        if isfield(l, 'read_noise_variance') ...
                && isfield(l, 'thermal_noise_variance_per_s') ...
                && isfield(l, 'pixel_baseline') ...
                && isfield(l, 'thermal_counts_per_s')
            [offsetmap, varmap] = makeExpDependMap(l, exposuretime_data);
        else
            offsetmap=l.offsetmap;
            varmap=l.varmap;
        end
        %assume gainmap is independent of exposure time
            gainmap=l.gainmap;

        %varmap has to be in units of photons^2
            varmap=varmap.*gainmap.^2; 
                %???? in units of photons, to be used later directly. 
                %A: Yes, conversion to photons^2 is not considered in GPUmleFIT
    else
        gainmap=[];offsetmap=[];varmap=[];
    end
    
end

function [offsetmap, varmap] = makeExpDependMap(l, exposuretime_data)
        exposure_time_s = exposuretime_data /1e3; % maps have units seconds, SMAP handles timing in ms
        
        varmap = l.read_noise_variance + l.thermal_noise_variance_per_s * exposure_time_s;
        offsetmap = l.pixel_baseline + l.thermal_counts_per_s * exposure_time_s;
end % function