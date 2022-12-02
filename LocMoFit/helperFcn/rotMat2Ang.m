function [xrot,yrot,zrot] = rotMat2Ang(R, par)
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