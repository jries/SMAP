function imout=transformImage(transformf,image,cam_pixnm,roi)
    if nargin<4
        roi=[0 0 size(image)];
    end
    sizeim=size(image);
    extxnm=[roi(1) roi(1)+roi(3)]*cam_pixnm(1)/1000;
    extynm=[roi(2) roi(2)+roi(4)]*cam_pixnm(end)/1000;
%             R = imref2d(sizeim,cam_pixnm,cam_pixnm);
    R = imref2d(sizeim,extxnm,extynm);
    imout = imwarp(image,R,transformf,'OutputView',R);
end