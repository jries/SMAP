% imP = im2pol(imC)
% ---------------------------------------
%
% Transform an image from carthesian to polar coordinate
%
% Inputs:
%  imC        	Input image in carthesian coordinate
%
% Outputs:
%  imP        	Output image in polar coordinate
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

function imP = im2pol(imC)

rMin=0; 
rMax=1;

[Ny,Nx] = size(imC); 
xc = (Ny+1)/2; yc = (Nx+1)/2; 
sx = (Ny-1)/2; sy = (Nx-1)/2;

Nr = 2*Ny; Nth = 2*Nx;

dr = (rMax - rMin)/(Nr-1); 
dth = 2*pi/Nth;


r = rMin:dr:rMin+(Nr-1)*dr; 
th = (0:dth:(Nth-1)*dth)'; 
[r,th]=meshgrid(r,th); 

x = r.*cos(th); y = r.*sin(th); 
xR = x*sx + xc; yR = y*sy + yc;

imP = interp2(imC, xR, yR,'cubic'); 

imP(isnan(imP)) = 0;

end