% rsc = linmap(val,valMin,valMax,mapMin,mapMax)
% ---------------------------------------
%
% Performs a linear mapping of val from the range [valMin,valMax] to the range [mapMin,mapMax]
%
% Inputs:
%  val        	Input value
%  valMin		Minimum value of the range of val
%  valMax		Maximum value of the range of val
%  mapMin		Minimum value of the new range of val
%  mapMax		Maximum value of the new range of val
%
% Outputs:
%  rsc        	Rescaled value
%
% Example : rsc = linmap(val,0,255,0,1); % map the uint8 val to the range [0,1]
%
% ---------------------------------------
%
% A detailled description of the method can be found in : 
% "Descloux, A., K. S. Grußmayer, and A. Radenovic. "Parameter-free image 
% resolution estimation based on decorrelation analysis."
% Nature methods (2019): 1-7."
% 
%   Copyright © 2018 Adrien Descloux - adrien.descloux@epfl.ch, 
%   Ecole Polytechnique Federale de Lausanne, LBEN,
%   BM 5.134, Station 17, 1015 Lausanne, Switzerland.
%
%  	This program is free software: you can redistribute it and/or modify
%  	it under the terms of the GNU General Public License as published by
% 	the Free Software Foundation, either version 3 of the License, or
%  	(at your option) any later version.
%
%  	This program is distributed in the hope that it will be useful,
%  	but WITHOUT ANY WARRANTY; without even the implied warranty of
%  	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  	GNU General Public License for more details.
%
% 	You should have received a copy of the GNU General Public License
%  	along with this program.  If not, see <http://www.gnu.org/licenses/>.

function  rsc = linmap(val,valMin,valMax,mapMin,mapMax)

if nargin == 3 % normalize the data between valMin and valMax
    mapMin = valMin;
    mapMax = valMax;
    valMin = min(val(:));
    valMax = max(val(:));
end

% convert the input value between 0 and 1
tempVal = (val-valMin)./(valMax-valMin);

% clamp the value between 0 and 1
map0 = tempVal < 0;
map1 = tempVal > 1;
tempVal(map0) = 0;
tempVal(map1) = 1;

% rescale and return
rsc = tempVal.*(mapMax-mapMin) + mapMin;

end