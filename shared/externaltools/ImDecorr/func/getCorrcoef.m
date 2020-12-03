% cc = getCorrcoef(I1,I2,c1,c2)
% ---------------------------------------
%
% Return the normalized correlation coefficient expressed in Fourier space
%
% Inputs:
%  I1        	Complex Fourier transfom of image 1
%  I2           Complex Fourier transfom of image 2
%
% Outputs:
%  cc        	Normalized cross-correlation coefficient
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

function cc = getCorrcoef(I1,I2,c1,c2)

if nargin < 4
    c2 = sqrt(sum(sum(abs(I2).^2)));
end
if nargin < 3
	c1 = sqrt(sum(sum(abs(I1).^2)));
end

cc = sum(sum(real(I1.*conj(I2))))./((c1*c2));
cc = floor(1000*cc)/1000;