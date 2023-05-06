function valid = spacewin_isvalid(W)
% SPACEWIN_ISVALID check whether a spatial window is valid
% L = SPACEWIN_ISVALID(W) check that the space window (ROI) is valid/allowable
%       W should be a struct with a 'type' field and others depending on
%       the type.
%       type = 'polygon'    A polygonal window, specified by fields 'x' and
%       'y' in the same units as the dataset that this window is to be used
%       for. x and y should be of the same size, and the ith vertex of the
%       polygonal mask is x(i), y(i).
%       type = 'polyshape'  A window of type polyshape, which is stored in
%       field 'p'
%       type = 'image'      A logical image specified by an imref2d in
%       field 'ref' and a 2d matrix in field 'im'. The ImageSize attribute
%       of 'ref' should match the size of 'im'. Matlab axis conventions
%       apply: the first index corresponds to the y axis, and the second
%       index corresponds to x.

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

valid = false;

if isstruct(W) && isfield(W,'type')
    switch W.type
        case 'polygon'
            valid = poly_isvalid(W);
        case 'polyshape'
            valid = polyshape_isvalid(W);
        case 'image'
            valid = image_isvalid(W);
        otherwise
            warning('spatial window type not known');
    end
end

end

function valid = poly_isvalid(W)
valid = isfield(W, 'x') && isfield(W,'y');
if valid
    valid = isequal(size(W.x), size(W.y)) && isvector(W.x);
    if ~valid
        warning('spacewin_isvalid: ''x'' and ''y'' should be vectors of the same size');
    end
else
    warning('spacewin_isvalid: spatial window with type ''polygon'' should have fields ''x'' and ''y''');
end
end

function valid = polyshape_isvalid(W)
valid = isfield(W, 'p') && isa(W.p, 'polyshape');
end

function valid = image_isvalid(W)
valid = isfield(W, 'im') && isfield(W, 'ref') && isa(W.im, 'logical') && isa(W.ref, 'imref2d');
if ~valid
    warning('spacewin_isvalid: spatial window with type ''im'' should have fields ''im'' and ''ref'' of type logical and imref2d, respectively');
end
end
