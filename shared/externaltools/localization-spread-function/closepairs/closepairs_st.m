% Copyright (C) 2021 Thomas Shaw, Frank Fazekas and Sarah Veatch
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

function [dxout,dyout,dtout,err] = closepairs_st_binned(x, y, t, rmax, taumin, taumax, noutmax)

ind = sortbywhich(x,y, rmax, t, taumin, taumax);
switch ind
    case 1
        % Sort based on x
        [x, s] = sort(x(:));
        
        t = t(s);
        y = y(s);
        sortbyt = uint64(0);
    case 2
        % Sort based on y (and swap x and y)
        tmp = y;
        y = x;
        x = tmp;
        [x, s] = sort(x(:));
        
        t = t(s);
        y = y(s);
        sortbyt = uint64(0);
    case 3
        % Sort based on t
        [t, s] = sort(t(:));
        
        x = x(s);
        y = y(s);
        sortbyt = uint64(1);
end

[dxout,dyout,dtout,err] = Fclosepairs_st(x,y,t, rmax, taumin, taumax, uint64(noutmax), sortbyt);
