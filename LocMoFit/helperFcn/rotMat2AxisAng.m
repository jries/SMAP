function k = rotMat2AxisAng(R)
% :func:rotMat2AxisAng calcualtes the corresponding rotation axis based on the rotation matrix R. 
%
% Uasage:
%   k = rotMat2AxisAng(R)
%
% Args:
%   R (numeric matrix): a rotation matrix.
%
% Returns:
%   k (numeric vector): a 1-by-3 vector. ||k|| is the rotation angle theta about the unit vector k/||k||.
%
% Last update:
%   06.12.2022
%
    theta = acos((trace(R)-1)./2); % this is also the norm of the vector k
    if theta == 0
        k = [0 0 0];
    elseif 0<theta&&theta<pi
        k = [R(3,2)-R(2,3) R(1,3)-R(3,1) R(2,1)-R(1,2)];
    elseif theta == pi
        [~,ind_max] = max([R(1,1) R(2,2) R(3,3)]);
        switch ind_max
            case 1
                k = [R(1,1)+1, R(1,2), R(1,3)];
            case 2
                k = [R(2,1), R(2,2)+1, R(2,3)];
            otherwise
                k = [R(3,1), R(3,2), R(3,3)+1];
        end
    end
    if theta ~= 0
        k = (k./norm(k)).*theta;
    end
end