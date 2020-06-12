% [im,pps] = smlmHist(data,pps,FOV)
% ---------------------------------------
%
% Render the localizations data as a bilinear histogram
%
% Inputs:
%  data        	Input localizations, should be of size [N,2]
%                   with N: Number of loc
%  pps          Target pixel size, should be in same units as data
%  FOV          All localizations within the range [0, FOV] will be rendered
%                   FOV should be in same units as data and pps
%
% Outputs:
%  im        	Rescaled value
%
% ---------------------------------------
%   Adrien Descloux - adrien.descloux@epfl.ch, 
%   Ecole Polytechnique Federale de Lausanne, LBEN,
%   BM 5.134, Station 17, 1015 Lausanne, Switzerland.

function [im,pps] = smlmHist(data,pps,FOV)
% input handling
if nargin < 3
FOVx = max(data(:,1));
FOVy = max(data(:,2));
else
    FOVx = FOV;
    if numel(FOV) == 1
        FOVy = FOV;
    else
        FOVy = FOV(2);
    end
end

% compute the number of pixels in the image
Nx = ceil(FOVx/pps);
Ny = ceil(FOVy/pps);
% if FOV is not a integer multiple of pps, adapt the pixel size
pps = FOVx/Nx;
FOVy = pps*Ny;

% map data in image space
data(:,1) = linmap(data(:,1),0,FOVx,1,Nx);
data(:,2) = linmap(data(:,2),0,FOVy,1,Ny);
% filter data out of the range [0, FOV]
x0 = floor(data(:,1)); y0 = floor(data(:,2));
map = x0 < Nx & x0 > 0 & y0 < Ny & y0 > 0;
x0 = x0(map); y0 = y0(map);
% compute linear weights
wx = data(map,1)-x0; wy = data(map,2)-y0;
wx = reshape(wx,[1 1 numel(wx)]);
wy = reshape(wy,[1 1 numel(wy)]);
% bilinear weights
M = bsxfun(@times,[1-wy; wy],[1-wx, wx]);
% rendering loop
im = zeros(Ny,Nx);
for k = 1:size(y0,1)
    im(y0(k),x0(k))     = im(y0(k),x0(k))       + M(1,1,k);
    im(y0(k),x0(k)+1)   = im(y0(k),x0(k)+1)     + M(1,2,k);
    im(y0(k)+1,x0(k))   = im(y0(k)+1,x0(k))     + M(2,1,k);
    im(y0(k)+1,x0(k)+1) = im(y0(k)+1,x0(k)+1)   + M(2,2,k);
end
