function [x, y, z, rotMat] = rotcoord3(x,y,z,xrot,yrot,zrot, orderRot)
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
    switch orderRot
        case 'XYZ'
            s1 = sin(xrot);
            s2 = sin(yrot);
            s3 = sin(zrot);
            c1 = cos(xrot);
            c2 = cos(yrot);
            c3 = cos(zrot);
            a11 = c2.*c3;
            a12 = -c2.*s3;
            a13 = s2;
            a21 = c1.*s3+c3.*s1.*s2;
            a22 = c1.*c3-s1.*s2.*s3;
            a23 = -c2.*s1;
            a31 = s1.*s3-c1.*c3.*s2;
            a32 = c3.*s1+c1.*s2.*s3;
            a33 = c1.*c2;
          case 'ZYX'
            s1 = sin(zrot);
            s2 = sin(yrot);
            s3 = sin(xrot);
            c1 = cos(zrot);
            c2 = cos(yrot);
            c3 = cos(xrot);
            a11 = c1.*c2;
            a12 = c1.*s2.*s3-c3.*s1;
            a13 = s1.*s3 + c1.*c3.*s2;
            a21 = c2.*s1;
            a22 = c1.*c3+s1.*s2.*s3;
            a23 = c3.*s1.*s2 - c1.*s3;
            a31 = -s2;
            a32 = c2.*s3;
            a33 = c2.*c3;
    end
    rotMat = zeros(4,4,length(a12));
    rotMat(1,1,:) = a11;
    rotMat(1,2,:) = a12;
    rotMat(1,3,:) = a13;
    rotMat(2,1,:) = a21;
    rotMat(2,2,:) = a22;
    rotMat(2,3,:) = a23;
    rotMat(3,1,:) = a31;
    rotMat(3,2,:) = a32;
    rotMat(3,3,:) = a33;
    rotMat(4,4,:) = ones(size(s1));
    locsTMat = pagemtimes(rotMat, locsMat);
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