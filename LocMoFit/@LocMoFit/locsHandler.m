function [locs, flag] = locsHandler(obj, locs, lParsVal, modelID, varargin)%!!! modelID needs to be remove.
%% LOCSHANDLER Handle the transformation of locs
% Transform the localizations based on the lPars
% This is applied to the localizations when the model is in the form of a
% image. This transformation equals the reverse operation acting on the model.
%
% Usage:
%   locsT = locsHandler(obj, locs, lParsVal)
% Arg:
%   locs: a struct of the localization table.
%   lParsVal: values of localization parameters.
%   locsT: a struct of the transformed locs.
% Note:
%	!!! locsHandler does not only handle locs anymore. The function name
%	should be changed to transformCoord.

%% Initiation
p = inputParser;
p.addParameter('usedformalism','rotationMatrix');
p.addParameter('target','locs');                % can be either 'locs' or 'model'; not used anymore
p.addParameter('order_transform','TR');
p.addParameter('onlyLocpre',false);

parse(p, varargin{:});
% need to implement scaling and weight            % !!! this might not be necessary
zrot = deg2rad(lParsVal.zrot);                      % from degree to radians
if obj.dataDim == 3
    xrot = deg2rad(lParsVal.xrot);                  % the same
    yrot = deg2rad(lParsVal.yrot);
end

if isfield(locs, 'locprecnm')
    locs.locprecnm = sqrt(locs.locprecnm.^2+lParsVal.variation.^2);
    if obj.dataDim == 3
        if any(locs.locprecznm==0)
            locs.locprecznm = locs.locprecnm*3;
            disp('locprecznm is missing: locprecnm*3 is used as locprecznm.')
        end
        locs.locprecznm = sqrt(locs.locprecznm.^2+lParsVal.variation.^2);
    end
    if p.Results.onlyLocpre
        flag = 1;
        return
    end
end

% shift the coordinates
x = locs.xnm;
y = locs.ynm;
if obj.dataDim == 3
    z = locs.znm;
end

switch p.Results.order_transform
    case 'TR'
        allActions = {'translation','scaling','rotation'};
        sign = -1;
    case 'RT'
        allActions = {'rotation','scaling','translation'};
        sign = 1;
end

%% Loop through twice: one for translation, on for rotation. The order of
% these two depends on the cellstr allActions.
l = 1;
while l <= length(allActions)
    action = allActions{l};
switch action
% Translation
    case 'translation'
        x = x+sign*lParsVal.x;
        y = y+sign*lParsVal.y;
        if obj.dataDim == 3
            z = z+sign*lParsVal.z;
        end
%Scaling
    case 'scaling'
        x = x.*(lParsVal.xscale.^sign);
        y = y.*(lParsVal.yscale.^sign);
        if obj.dataDim == 3
            z = z.*(lParsVal.zscale.^sign);
        end
%Rotation
    case 'rotation'
        usedformalism = p.Results.usedformalism;
        if obj.dataDim == 2
            [x,y] = rotcoord2(x,y,sign*zrot');
        else
            switch usedformalism
                case 'legacy'
                    if obj.dataDim == 3
                        xrot = deg2rad(lParsVal.xrot);                  % the same
                        yrot = deg2rad(lParsVal.yrot);
                        z = locs.znm;
                        [x,z] = rotcoord2(x,z,-yrot');
                        [y,z] = rotcoord2(y,z,-xrot');
                    end
                    flag = 1;
                case 'rotationMatrix'
                    % extrinsic rotation (ZYX) (order: X->Y->Z)
                    [x,y,z] = rotcoord3(x, y, z, sign*xrot, sign*yrot, sign*zrot, 'ZYX');
                    x = x';
                    y = y';
                    z = z';
                    flag = 1;
                case 'rotationMatrixRev'
                    % extrinsic rotation (XYZ) (order: Z->Y->X)
                    % this one is used when the locsHandler is applied to model
                    % coordinates.
                    [x,y,z] = rotcoord3(x, y, z, sign*xrot, sign*yrot, sign*zrot, 'XYZ');
                    x = x';
                    y = y';
                    z = z';
                    flag = 1;
                case 'quaternion'
                    qk = lParsVal.zrot;
                    qi = lParsVal.xrot; qj = lParsVal.yrot;
                    du = sum([qi qj qk].^2);
                    if du > 1
                        flag = 0;
                    else
                        qw = (double(1)-du).^0.5;
                        rotMat = quaternionRotMat(qw, qi, qj, qk);
                        rotAng = quaternion2deg(qw,qi, qj, qk);

                        rot = {'xrot','yrot','zrot'};
                        indFinal = true;

                        k = 0;
                        % internal validation based on the rotational angles
                        while indFinal == true&&k<3
                            k = k+1;
                            [~,ind] = obj.getVariable(['par.m1.lPar.' rot{k}]);
                            if obj.allParsArg.fix(ind)
                                tol = 1;% tolerance in degree
                                if rotAng(k) > obj.allParsArg.value(ind)-tol && rotAng(k) < obj.allParsArg.value(ind)+tol
                                    indFinal = indFinal&&true;
                                else
                                    indFinal = false;
                                end
                            else
                                % follow the same rules as for the other parameters
                                if obj.allParsArg.value(ind)+obj.allParsArg.lb(ind) < obj.allParsArg.min(ind)
                                    minVal = obj.allParsArg.min(ind);
                                else
                                    minVal = obj.allParsArg.value(ind)+obj.allParsArg.lb(ind);
                                end
                                if obj.allParsArg.value(ind)+obj.allParsArg.ub(ind) > obj.allParsArg.max(ind)
                                    maxVal = obj.allParsArg.max(ind);
                                else
                                    maxVal = obj.allParsArg.value(ind)+obj.allParsArg.ub(ind);
                                end

                                if ~(rotAng(k)>minVal && rotAng(k)<maxVal)
                                    indFinal = false;
                                end
                            end
                        end
                        if indFinal
                            locsMat = zeros(3,length(x));
                            locsMat(1,:) = x;
                            locsMat(2,:) = y;
                            locsMat(3,:) = z;
                            locsTMat = ndfun('mult',rotMat, locsMat);
                            x = locsTMat(1,:)';
                            y = locsTMat(2,:)';
                            z = locsTMat(3,:)';
                            flag = 1;
                        else
                            flag = 0;
                        end
                    end
            end
        end
end
l = l+1;
end
%% Output
if obj.dataDim == 3
    locs.znm = z;
end
locs.xnm = x;
locs.ynm = y;
%% Variation
if isfield(locs,'locprecnm')&&any(lParsVal.variation>0)
    locs.locprecnm = sqrt(locs.locprecnm.^2+lParsVal.variation.^2);
    if isfield(locs,'locprecznm')
        locs.locprecznm = sqrt(locs.locprecznm.^2+lParsVal.variation.^2);
    end
end
end

%% Internal functions
% Quaternion-derived rotation matrix
function rotMat = quaternionRotMat(qw,qi,qj,qk)
    s = 1;
    m11 = 1-2.*s.*(qj.^2+qk.^2);
    m12 = 2.*s.*(qi.*qk-qk.*qw);
    m13 = 2.*s.*(qi.*qk+qj.*qw);
    m21 = 2.*s.*(qi.*qk+qk.*qw);
    m22 = 1-2.*s.*(qi.^2+qk.^2);
    m23 = 2.*s.*(qi.*qk-qi.*qw);
    m31 = 2.*s.*(qi.*qk-qj.*qw);
    m32 = 2.*s.*(qj.*qk+qi.*qw);
    m33 = 1-2.*s.*(qi.^2+qj.^2);
    rotMat = zeros(3,3,length(qw));
    rotMat(1,1,:) = m11;
    rotMat(1,2,:) = m12;
    rotMat(1,3,:) = m13;
    rotMat(2,1,:) = m21;
    rotMat(2,2,:) = m22;
    rotMat(2,3,:) = m23;
    rotMat(3,1,:) = m31;
    rotMat(3,2,:) = m32;
    rotMat(3,3,:) = m33;
end