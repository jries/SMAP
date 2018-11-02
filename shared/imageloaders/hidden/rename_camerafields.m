file='settings/cameras.mat'; 
load(file)
for k=1:length(cameras)
    for l=1:length(cameras(k).state)
        cameras(k).state(l).par{2,1}='cam_pixelsize_um';
    end;
end
for k=1:length(cameras)
    cameras(k).par{2,1}='cam_pixelsize_um';
end

save(file,'cameras','camtab')