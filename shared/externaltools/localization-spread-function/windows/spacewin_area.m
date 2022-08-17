function A = spacewin_area(W)
% SPACEWIN_AREA compute area of a spatial window
% A = SPACEWIN_AREA(W)      W should be a valid spatial window according to
%                           SPACEWIN_ISVALID(). A is the area of the
%                           spatial window, in appropriate units.

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

if ~spacewin_isvalid(W)
    error('spacewin_area: invalid spatial window')
end

switch W.type
    case 'polygon'
        A = polyarea(W.x, W.y);
    case 'polyshape'
        A = W.p.area();
    case 'image'
        A = sum(W.im(:))*W.ref.PixelExtentInWorldX*W.ref.PixelExtentInWorldY;
    otherwise
        error('spacewin_area: not implemented for this type of spatial window');
end
