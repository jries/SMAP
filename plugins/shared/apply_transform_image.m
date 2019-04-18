function imout=apply_transform_image(img,transform,p)
if isa(transform,'interfaces.LocTransformN')
    switch p.datapart.selection
        case {'reference','all (R->T)'}
        imout=transform.transformImageToTarget(2,img,p.cam_pixelsize_nm,p.roitiff);
        case {'all (T->R)','target'}
        imout=transform.transformImageToReference(2,img,p.cam_pixelsize_nm,p.roitiff);
        otherwise
        disp('wrong transformation specified')  
    end
else
switch p.datapart.selection
    case {'reference','all (R->T)'}
    imout=transform.transformImageFwd(img,p.cam_pixelsize_nm,p.roitiff);
    case {'all (T->R)','target'}
    imout=transform.transformImageInv(img,p.cam_pixelsize_nm,p.roitiff);
    otherwise
        disp('wrong transformation specified')      
end
end

