function R = rodringues2rotMat(k, theta)
% :func:`rodringues2rotMat` calcualtes the corresponding rotation matrix :attr:`R` based on the rotation-axis :attr:`k` and angle :attr:`theta`.
%
% Uasage:
%   R = rodringues2rotMat(k, theta)
%
% Args:
%   k (numeric vector): a 1-by-3 unit vector. This is the rotation axis vector with ||k|| = 1.
% 	theta (numeric scalar): rotation angle about the unit vector k.
%
% Returns:
%   R (numeric matrix): a 3-by-3 rotation matrix.
%
% Last update:
%   06.12.2022
%
    S = [0 -k(3) k(2);...
        k(3) 0 -k(1);...
        -k(2) k(1) 0];
    I = diag(ones([1,3]));
%     R = (I+sin(theta).*S+(1-cos(theta).*S.^2));
    R = cos(theta).*I+(1-cos(theta)).*k.*k'+sin(theta).*S;
end


