% [out,mask] = apodImRect(in,N)
% ---------------------------------------
%
% Apodize the edges of a 2D image
%
% Inputs:
%  in        	Input image
%  N            Number of pixels of the apodization
%
% Outputs:
%  out        	Apodized image
%  mask         Mask used to apodize the image
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

function [out,mask] = apodImRect(in,N)

Nx = size(in,2);
Ny = size(in,1);

x = abs(linspace(-Nx/2,Nx/2,Nx));
y = abs(linspace(-Ny/2,Ny/2,Ny));
mapx = x > Nx/2 - N;
mapy = y > Ny/2 - N;

val = mean(mean(in(:)));

d = (-abs(x)- mean(-abs(x(mapx)))).*mapx;
d = linmap(d,-pi/2,pi/2);
d(not(mapx)) = pi/2;
maskx = (sin(d)+1)/2;

d = (-abs(y)- mean(-abs(y(mapy)))).*mapy;
d = linmap(d,-pi/2,pi/2);
d(not(mapy)) = pi/2;
masky = (sin(d')+1)/2;

% make it 2D
mask = repmat(masky,[1 size(in,2)]).*...
    repmat(maskx,[size(in,1) 1]);

out = (in-val).*mask + val;