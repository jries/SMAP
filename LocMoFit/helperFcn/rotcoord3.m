function [x, y, z, R] = rotcoord3(x,y,z,xrot,yrot,zrot, orderRot)
% The orderRot is the order of rotation.
% We are using the Tait–Bryan angles here.
% The angle here represents an extrinsic rotation.
lMultiDim = ~all([size(x,2) size(y,2) size(z,2) size(xrot,2) size(yrot,2) size(zrot,2)]==1);
    space = zeros(size(x));
    if ~lMultiDim
        locsMat = zeros(4,size(x,1));
        locsMat(1,:) = x;
        locsMat(2,:) = y;
        locsMat(3,:) = z;
        locsMat(4,:) = ones(size(x));
    else
        locsMat = zeros(4,size(x,1),size(x,2));
        locsMat(1,:,:) = x+space;
        locsMat(2,:,:) = y+space;
        locsMat(3,:,:) = z+space;
        locsMat(4,:,:) = ones(size(x));
    end
    R = rotAng2rotMat(xrot, yrot, zrot, orderRot);
    locsTMat = pagemtimes(R, locsMat);
%     locsTMat = ndfun('mult',rotMat, locsMat);
    if ~lMultiDim
        x = locsTMat(1,:);
        y = locsTMat(2,:);
        z = locsTMat(3,:);
    else
        x = squeeze(locsTMat(1,:,:))';
        y = squeeze(locsTMat(2,:,:))';
        z = squeeze(locsTMat(3,:,:))';
    end
end