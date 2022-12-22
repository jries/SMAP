function rotMat = rotAng2rotMat(xrot,yrot,zrot, orderRot)
% Get the rotation matrix based on Tait–Bryan angles.
% :func:`rotAng2rotMat` calcualtes the corresponding rotation matrix :attr:`rotMat` based on Tait–Bryan angles, given a specific order.
%
% Uasage:
%   rotMat = rotAng2rotMat(xrot,yrot,zrot, orderRot)
%
% Args:
%   xrot (numeric scalar): rotation angle about the x-axis.
%	yrot (numeric scalar): rotation angle about the y-axis.
%	zrot (numeric scalar): rotation angle about the z-axis.
%	orderRot (string): rotation order. Either 'XYZ' or 'ZYX'. It is the order of the respective matrix derived from each rotation angle.
%
% Returns:
%   rotMat (numeric matrix): a 4-by-4 translation matrix.
%
% Last update:
%   06.12.2022
%
switch orderRot
    case 'XYZ'
        % extrinsic rotation (XYZ) (order: X->Y->Z)
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
        % xtrinsic rotation (XYZ) (order: Z->Y->X)
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
end