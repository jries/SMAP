function [xrot,yrot,zrot] = rotMat2Ang(R, par)
% :func:`rotMat2Ang` calcualtes the corresponding Tait–Bryan angles based on the rotation matrix :attr:`R`, given a specific parameterization :attr:`par`.
%
% Uasage:
%   [xrot,yrot,zrot] = rotMat2Ang(R, par)
%
% Args:
%   R (numeric matrix): a 3-by-3 rotation matrix.
% 	par (string): either 'rotationMatrix' or 'rotationMatrixRev'.
%
% Returns:
%   xrot (numeric scalar): rotation angle about the x-axis.
%	yrot (numeric scalar): rotation angle about the y-axis.
%	zrot (numeric scalar): rotation angle about the z-axis.
%
% Last update:
%   06.12.2022
%
    switch par
        case 'rotationMatrix'
            % extrinsic rotation (ZYX) (order: X->Y->Z)
            switch R(1,3)
                case -1
                    yrot = pi/2;
                    xrot = atan2(R(2,1),R(3,1));
                    zrot = 0;
                case 1
                    yrot = -pi/2;
                    xrot = atan2(-R(2,1),-R(3,1));
                    zrot = 0;
                otherwise
                    yrot = asin(-R(1,3));
                    if -pi/2<yrot && yrot<pi/2
                        xrot = atan2(R(2,3),R(3,3));
                        zrot = atan2(R(1,2),R(1,1));
                    else
                        xrot = atan2(-R(2,3),-R(3,3));
                        zrot = atan2(-R(1,2),-R(1,1));
                    end
            end
        case 'rotationMatrixRev'
            % extrinsic rotation (XYZ) (order: Z->Y->X)
            switch R(3,1)
                case -1
                    yrot = pi/2;
                    xrot = atan2(R(1,2),R(1,3));
                    zrot = 0;
                case 1
                    yrot = -pi/2;
                    xrot = atan2(-R(1,2),-R(1,3));
                    zrot = 0;
                otherwise
                    yrot = asin(-R(3,1));
                    if -pi/2<yrot && yrot<pi/2
                        xrot = atan2(R(3,2),R(3,3));
                        zrot = atan2(R(2,1),R(1,1));
                    else
                        xrot = atan2(-R(3,2),-R(3,3));
                        zrot = atan2(-R(2,1),-R(1,1));
                    end
            end
    end
    xrot = rad2deg(xrot);
    yrot = rad2deg(yrot);
    zrot = rad2deg(zrot);
end