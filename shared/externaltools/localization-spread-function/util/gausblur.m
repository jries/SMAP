function I2 = gausblur(I1, sigma, cutoff)
% B = GAUSBLUR(A, SIGMA)
%
% A     Image (2d double array) to be blurred
% sigma Sigma of gaussian to blur by, in image pixel units
%
% B     Output image, with same size as A

% Copyright (C) 2021 Thomas Shaw
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
    % this is about 99th percentile.
    cutoff = 3;
end

mxx = ceil(cutoff*sigma);
sz = 2*mxx + 1;

% Pixelate by averaging the gaussian over the pixel:
h = diff(erf( ((-mxx-.5):(mxx+.5))/(sqrt(2)*sigma) ));
h = h/sum(h);

% Pixelate by evaluating the gaussian at the center of the pixel
%h = exp(-((-mxx):mxx).^2/2/sigma^2);
%h = h/sum(h);

I2 = conv2(h,h,I1,'same');
%I2 = conv2(I1, fspecial('gaussian', sz, sigma), 'same');
