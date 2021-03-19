function [probability, v1sum, v2sum, ampC, nV1, nV2] = nucleationNCoat3(x, y, z, xcenter, ycenter, distance, ampC, angle, innerRadR, outerRadR, innerRadC, outerRadC, thickness, varargin)
 
%% Rotate coordinates
[x,y] = rotcoord(x,y,-angle*pi/180);

%% Calculate volumes of structures
if outerRadR == 70 && innerRadR==30 && thickness==60            % Pre-calculated for common parameters
    volR = 1.5080e+06;
else                                                            % Real-time for uncommon parameters
    volR = (2*pi.*thickness).*(outerRadR.^2 - innerRadR.^2);
end

if outerRadC == 40&&innerRadC==20                               % Pre-calculated for common parameters
    volC = 1.1728e+05;
else                                                            % Real-time for uncommon parameters
    volC = (pi*2/3).*(outerRadC.^3 - innerRadC.^3);
end

%% Calculate the sum of density
% for normalizating the sum of density to one. The sum of the cutoff needs
% to be counted in. Therefore the roi size and the unit of offset have to
% be known.
if length(varargin) > 1                                         % Checking for name&value pairs
    indROISize = find(arrayfun(@(x)isequal({'RoiSize'}, x), varargin));
    indOffset = find(arrayfun(@(x)isequal({'SumOffset'}, x), varargin));
    if indROISize > 0                                           % get roiSize
        roiSize = varargin{indROISize+1};
    end
    if indOffset > 0                                            % get offset
        SumOffset = varargin{indOffset+1};
        sumAllVol = 1+SumOffset;
        offset = SumOffset/(roiSize^2);
    else
        offset = 1.1111e-06;
        sumAllVol = 1.1;
    end
end

%% the main part
if length(z)>1                                                  % dual-colour mode
    ch1 = z==1;
    ch2 = ~ch1;
    xch1 = x(ch1); ych1 = y(ch1);
    xch2 = x(ch2); ych2 = y(ch2);
    v1 = 0.5.*thickRing([xcenter, ycenter, innerRadR, outerRadR, thickness], xch1, ych1)./volR;
    v2 = 0.5.*cap([xcenter, ycenter+distance, innerRadC, outerRadC, thickness], xch2, ych2)./volC;
    ampR = 1;
    if ampC == 0
        nV1 = sum(v1 ~= 0);
        nV2 = sum(v2 ~= 0);
        if(nV1)>0
            ampC = nV2/nV1;
        else
            ampC = offset;
        end
    end
    ampTotal = ampR+ampC;
    v1 = v1.*ampR./ampTotal;
    v2 = v2.*ampC./ampTotal;
    v1sum = sum(v1); v2sum = sum(v2);
%     totalValue = v1sum+v2sum;
%     v1 = v1.*v1sum./totalValue;
%     v2 = v2.*v2sum./totalValue;
    v =[v1;v2];
else                                                            % single-colour mode
    % 0.5 * total volue given a pair of x and y / total volume
    % = 0.5 * desity given a pair of x and y
    % the overall density of the unit structure will be 1
    v1 = 0.5.*thickRing([xcenter, ycenter, innerRadR, outerRadR, thickness], x, y)./volR;
    v2 = 0.5.*cap([xcenter, ycenter+distance, innerRadC, outerRadC, thickness], x, y)./volC;
    ampR = 1;
    if ampC == 0
        nV1 = sum(v1 ~= 0);
        nV2 = sum(v2 ~= 0);
        if(nV1)>0
            ampC = nV2/nV1;
        else
            ampC = offset;
        end
    end
    nV1 = sum(v1 ~= 0);
    nV2 = sum(v2 ~= 0);
    ampTotal = ampR+ampC;
    v1 = v1.*ampR./ampTotal;
    v2 = v2.*ampC./ampTotal;
    v1sum = sum(v1); v2sum = sum(v2);
    v = v1+v2;
end
probability = (v + offset)./sumAllVol;
end