function correction = spatial_edge_correction(W, r, hy)
% SPATIAL_EDGE_CORRECTION spatial edge correction for a correlation function
% of data observed on a given spatial window
% C = SPATIAL_EDGE_CORRECTION(W,R)      Isotropic spatial edge correction
%             for pair correlation functions or Ripley's K functions, for a
%             spatial window W (as specified for spacewin_isvalid()), at
%             distances R
%             Normalized to Area, so that it is order 1.
% C = SPATIAL_EDGE_CORRECTION(W,HX,HY)  Anisotropic spatial edge correction.
%             As above, but for a vector displacement hx,hy in x and y
%             respectively.

% Copyright (C) 2021 Thomas Shaw, and Sarah Veatch
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

if nargin < 3
    isotropic = true;
else
    isotropic = false;
    hx = r;
    if ~isequal(size(hx), size(hy))
        error('spatial_edge_correction: for anisotropic edge correction, hx and hy must be same size')
    end
end

if strcmp(W.type, 'image')
    % Strategy is a bit different for images, for computational reasons
    wij = conv2_fft_norm_and_interp(W.im, W.ref);
else
    area = spacewin_area(W);
    wij = @(dx,dy) wij_poly(W,dx,dy)/area;
end

if isotropic
    ntheta = 60;
    thetas = pi*(1:ntheta)/ntheta; % only need 0 to pi because w_ij(theta + pi) = w_ij(theta)
    correction = zeros(size(r));
    for i = 1:numel(r)
        for j = 1:ntheta
            dx = cos(thetas(j))*r(i);
            dy = sin(thetas(j))*r(i);
            correction(i) = correction(i) + wij(dx, dy)/ntheta;
        end
    end
    correction = correction(:);
else %anisotropic
    correction = zeros(size(hx));
    for i=1:numel(hx) % possibly this could be vectorized for the image case. I don't think it can be for the polygon case
        correction(i) = wij(hx(i),hy(i));
    end
end
end

function w = wij_poly(W, dx, dy)
switch W.type
    case 'polygon'
        mx = W.x;
        my = W.y;
        [x, y] = polybool('intersection', mx, my, mx + dx, my + dy);
        w = polyarea(x,y);
    case 'polyshape'
        p = W.p;
        pt = p.translate(dx, dy);
        pinter = intersect(p,pt);
        w = area(pinter);
end
end

function W = conv2_fft_norm_and_interp(I, ref)

sza = size(I, 1)*2 + 1;
szb = size(I, 2)*2 + 1;

C2 = abs(fftshift(ifft2(abs(fft2(I, sza, szb)).^2)))./sum(I(:));

x =  (-size(I,2):size(I,2))*ref.PixelExtentInWorldX;
y =  (-size(I,1):size(I,1))*ref.PixelExtentInWorldY;

[X,Y] = ndgrid(x,y);

W = griddedInterpolant(X,Y,C2');
end
