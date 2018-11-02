function loc=concentratefilelist(loc)
%file
for k=1:length(loc.file)
    try
    loc.file.info.allmetadata.files=char(loc.file.info.allmetadata.files');
    catch
         disp('files')
    end
end
%history
for k=1:length(loc.history)
    try
        loc.history{k}.children.CameraConverter.loc_cameraSettings.allmetadata.files=char(loc.history{k}.children.CameraConverter.loc_cameraSettings.allmetadata.files');
    catch
        
        try 
            loc.history(k).children.CameraConverter.loc_cameraSettings.allmetadata.files=char(loc.history(k).children.CameraConverter.loc_cameraSettings.allmetadata.files');
        catch
            disp('history')
        end
    end
end

for k=1:length(loc.siteexplorer.files)
    try
        loc.siteexplorer.files(k).info.allmetadata.files=char(loc.siteexplorer.files(k).info.allmetadata.files');
    catch
         disp('se')
    end
end

%siteexplorer
end


