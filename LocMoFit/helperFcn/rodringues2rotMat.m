%% Rodrigues% rotation formula

function R = rodringues2rotMat(k, theta)
    S = [0 -k(3) k(2);...
        k(3) 0 -k(1);...
        -k(2) k(1) 0];
    I = ones([1,3]);
    R = (I+sin(theta)*S+(1-cos(theta)*S^2));
end


