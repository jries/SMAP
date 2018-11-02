function imbwout=imdilate3(bwim,dilation)
se = strel('disk',dilation);
imbw3 = imdilate(bwim,se);
imbw3 = imdilate(permute(imbw3,[2 3 1]),se);
imbw3 = imdilate(permute(imbw3,[2 3 1]),se);
imbwout=permute(imbw3,[2 3 1]);

