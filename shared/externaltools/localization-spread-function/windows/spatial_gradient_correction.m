function [gm, gamma] = spatial_gradient_correction(x1,y1,x2,y2,spacewin,r, sigma, varargin)
% SPATIAL_GRADIENT_CORRECTION Compute correction factor to correct
% spatial cross correlations for spatially varying density.
% 
% [GM, GAMMA] = SPATIAL_GRADIENT_CORRECTION(X1,Y1,X2,Y2,SPACEWIN,R,SIGMA)
%               Compute correction for a cross-correlation of the datasets
%               x1,y1 and x2,y2, on a region of interest defined by
%               SPACEWIN, using a kernel density estimator for the density
%               with Gaussian blurring with standard deviation of SIGMA. R
%               sets the spatial separations at which to evaluate the
%               corrections. GAMMA is a 2d array with the (anisotropic)
%               correction factors as a function of dx, dy. GM is averaged
%               over angles so that GM(i) is the correction for the
%               (isotropic) cross-correlation at R(i).
%
% [_] = SPATIAL_GRADIENT_CORRECTION(_,'type', t)
%               As above, specifying the type of edge correction to perform
%               when calculating the correction. t must be either 'points'
%               (the default), or 'pairs'. The 'points' edge correction
%               computes the correction by convolving an edge corrected
%               kernel density estimate for each channel independently. This
%               is the method described in Shaw, Moller, Waagepetersen. 2020.
%               The 'pairs' edge correction computes the correction on the
%               blurred distribution of between-channel pair displacements.

% Copyright (C) 2022 Thomas Shaw, and Sarah Veatch
% This file is part of SMLM SPACETIME RESOLUTION
% SMLM SPACETIME RESOLUTION is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% SMLM SPACETIME RESOLUTION is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% You should have received a copy of the GNU General Public License
% along with SMLM SPACETIME RESOLUTION.  If not, see <https://www.gnu.org/licenses/>
p = inputParser;

addParameter(p, 'type', 'points',@(x) any(validatestring(x,{'points', 'pairs'})))
parse(p, varargin{:})
type = p.Results.type;

d = r(2)-r(1);

switch spacewin.type
    case 'polygon'
        left = min(spacewin.x);
        right = max(spacewin.x);
        top = min(spacewin.y);
        bottom = max(spacewin.y);
    case 'image'
        left = spacewin.ref.XWorldLimits(1);
        right = spacewin.ref.XWorldLimits(2);
        top = spacewin.ref.YWorldLimits(1);
        bottom = spacewin.ref.YWorldLimits(2);
    case 'polyshape'
        points = spacewin.p.Vertices;
        left = min(points(:, 1));
        right = max(points(:, 1));
        top = min(points(:, 2));
        bottom = max(points(:, 2));
end
I1 = histcounts2(x1, y1, left:d:right, top:d:bottom)';
I2 = histcounts2(x2, y2, left:d:right, top:d:bottom)';


[X, Y] = meshgrid(left+d/2:d:right-d/2, top+d/2:d:bottom-d/2);
mask = spacewin_isinside(X(:),Y(:),spacewin);
mask = reshape(mask, size(X));

lx = 2*size(I1, 1)+1;
ly = 2*size(I1, 2)+1;

switch type
    case 'points'
        Iblur1 =gausblur(I1.*mask, sigma/d).*mask;
        Iblur2 =gausblur(I2.*mask, sigma/d).*mask;
        Mblur =gausblur(mask, sigma/d).*mask;
        
        Iblur1 = Iblur1./Mblur;
        Iblur1(~Mblur) = 0;
        Iblur2 = Iblur2./Mblur;
        Iblur2(~Mblur) = 0;
        
        I1_fudgefact = mean(Iblur1(logical(mask(:))))/mean(I1(logical(mask(:))));
        I1_densfact = mean(Iblur1(logical(mask(:))));
        I2_fudgefact = mean(Iblur2(logical(mask(:))))/mean(I2(logical(mask(:))));
        I2_densfact = mean(Iblur2(logical(mask(:))));
        
        gamma = fftshift(ifft2(fft2(Iblur1, lx, ly).*conj(fft2(Iblur2, lx, ly))))./ ...
            fftshift(ifft2(abs(fft2(mask, lx, ly)).^2))/...
            I2_densfact/I1_densfact*I1_fudgefact*I2_fudgefact;
    
    case 'pairs'
        xc = ifft2(fft2(I1, lx, ly).*conj(fft2(I2, lx, ly)));
        if numel(data)==1
            xc(1,1) = xc(1,1) - sum(I1(:));
        end
        xc = fftshift(xc);
        
        m2 = fftshift(ifft2(abs(fft2(mask, lx, ly)).^2));
        m2(m2 < 1) = 0; % prevent some problems with rounding errors later.
        mask2 = double(logical(m2));
        edgem2 = gausblur(mask2, sqrt(2)*sigma/d).*mask2;
        
        rat = xc ./ m2;
        rat(~mask2) = 0;
        ratblur = gausblur(rat, sqrt(2)*sigma/d);
        
        I1_densfact = mean(I1(logical(mask(:))));
        I2_densfact = mean(I2(logical(mask(:))));
        
        
        gamma = ratblur./edgem2/I2_densfact/I1_densfact;
        gamma(~mask2) = 0;
    otherwise
        error('spatial_gradient_correction: unknown edge correction type ''%s''', type);
end


% angular average
[X, Y] = meshgrid((1:size(gamma, 2))-size(gamma, 2)/2-1, (1:size(gamma, 1))-floor(size(gamma, 1)/2)-1);
R = d*sqrt(X.^2+Y.^2);
gm = zeros(size(r));
for i=1:numel(r)
    inds = R(:)>=r(i)-d/2 &R(:)<r(i)+d/2;
    gm(i) = mean(gamma(inds));
end
gm = gm';
