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

function [counts] = crosspairs_st_binned(x1, y1, t1, x2, y2, t2, ...
        rmax, nrout, taumin, taumax, ntout)

ind = sortbywhich(x1,y1,rmax,t1,taumin,taumax);
switch ind
    case 1
        % Sort based on x
        [x1, s1] = sort(x1(:));
        [x2, s2] = sort(x2(:));

        t1 = t1(s1);
        y1 = y1(s1);
        t2 = t2(s2);
        y2 = y2(s2);

        sortbyt = uint64(0);
    case 2
        % Sort based on y, then swap x and y
        tmpx1 = x1;
        tmpx2 = x2;
        [x1, s1] = sort(y1(:));
        [x2, s2] = sort(y2(:));

        t1 = t1(s1);
        y1 = tmpx1(s1);
        t2 = t2(s2);
        y2 = tmpx2(s2);

        sortbyt = uint64(0);
    case 3
        % Sort based on t
        [t1, s1] = sort(t1(:));
        [t2, s2] = sort(t2(:));

        y1 = y1(s1);
        x1 = x1(s1);
        y2 = y2(s2);
        x2 = x2(s2);

        sortbyt = uint64(1);
end


[counts] = Fcrosspairs_st_binned(x1,y1,t1, x2, y2, t2,...
                                    rmax, uint64(nrout), taumin, taumax, uint64(ntout),sortbyt);

counts = double(counts);
