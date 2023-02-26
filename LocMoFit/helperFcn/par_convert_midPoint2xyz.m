function [ctrlX,ctrlY,ctrlZ] = par_convert_midPoint2xyz(par, numOfCtrlPointSet)
% control points:
    xMid = par.xMid;
    yMid = par.yMid;
    zMid = par.zMid;
    dist = par.dist;
    
    % left side
    for k = 1:numOfCtrlPointSet
        rotAziL(k,:) = par.(['rotAziL' num2str(k)]); % ->
        rotEleL(k,:) = par.(['rotEleL' num2str(k)]); % ->
        rotAziR(k,:) = par.(['rotAziR' num2str(k)]); % ->
        rotEleR(k,:) = par.(['rotEleR' num2str(k)]); % ->
    end
    rotAziL = cumsum(rotAziL,1);
    rotEleL = cumsum(rotEleL,1);
    rotAziL = deg2rad(rotAziL);
    rotEleL = deg2rad(rotEleL);
    [diffXL,diffYL,diffZL] = sph2cart(rotAziL,rotEleL,-dist);
    
    rotAziR = cumsum(rotAziR,1);
    rotEleR = cumsum(rotEleR,1);
    rotAziR = deg2rad(rotAziR);
    rotEleR = deg2rad(rotEleR);
    [diffXR,diffYR,diffZR] = sph2cart(rotAziR,rotEleR,dist);
    
    %% 191128
    diffXL = cumsum([xMid; diffXL],1);
    diffYL = cumsum([yMid; diffYL],1);
    diffZL = cumsum([zMid; diffZL],1);
    
    diffXR = cumsum([xMid; diffXR],1);
    diffYR = cumsum([yMid; diffYR],1);
    diffZR = cumsum([zMid; diffZR],1);
    
    ctrlX = [diffXL(end:-1:2);diffXR];
    ctrlY = [diffYL(end:-1:2);diffYR];
    ctrlZ = [diffZL(end:-1:2);diffZR];
end