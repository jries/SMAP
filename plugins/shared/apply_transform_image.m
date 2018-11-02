function imout=apply_transform_image(img,transform,p)
switch p.datapart.selection
    case {'reference','all (R->T)'}
    imout=transform.transformImageFwd(img,p.cam_pixelsize_nm,p.roitiff);
    case {'all (T->R)','target'}
    imout=transform.transformImageInv(img,p.cam_pixelsize_nm,p.roitiff);
    otherwise
        disp('wrong transformation specified')      
end
