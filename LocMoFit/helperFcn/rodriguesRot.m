function xyz = rodriguesRot(xyz, k)
    theta = norm(k);
    k = k./theta;
    nLocs = size(xyz,1);
    k = repmat(k,[nLocs 1]);

    xyz = xyz.*cos(theta)+cross(k,xyz).*sin(theta)+k.*(dot(k,xyz,2)).*(1-cos(theta));
end