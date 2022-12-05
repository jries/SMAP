%% Rodrigues% rotation formula

function R = rodringues2rotMat(k, theta)
    S = [0 -k(3) k(2);...
        k(3) 0 -k(1);...
        -k(2) k(1) 0];
    I = diag(ones([1,3]));
%     R = (I+sin(theta).*S+(1-cos(theta).*S.^2));
    R = cos(theta).*I+(1-cos(theta)).*k.*k'+sin(theta).*S;
end


