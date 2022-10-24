function ind = spacewin_isinside(x,y,W)
% SPACEWIN_ISINSIDE check which points are inside a spatial window
% IND = SPACEWIN_ISINSIDE(X,Y,W)
%       IND(i) is true if the point X(i), Y(i) is in the spatial window W,
%       as specified by spacewin_isvalid(), and false otherwise.

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
    error('spacewin_isinside: invalid spatial window provided');
end

switch W.type
    case 'polygon'
        ind = inpolygon(x,y,W.x, W.y);
    case 'polyshape'
        [interior, onedge] = isinterior(W.p, x,y);
        ind = interior | onedge;
        ind = reshape(ind,size(x));
    case 'image'
        [r,c] = W.ref.worldToSubscript(x,y);
        i1 = ~isnan(r) & ~isnan(c);
        j = sub2ind(size(W.im), r(i1),c(i1));
        i2 = logical(W.im(j));
        ind = i1;
        ind(i1) = i2;
    otherwise
        error('spacewin_inside: spatial window of this type is not yet supported')
end
