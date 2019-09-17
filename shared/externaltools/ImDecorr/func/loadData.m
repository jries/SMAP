% [im,path] = loadData(path)
% ---------------------------------------
%
% Generic function to load arbitrary image stacks in the workspace
%
% Inputs:
%  path        	Full path of the file, if not provided the function will call uigetfile
%
% Outputs:
%  im        	Loaded image
%  path 		Image full path
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

function [im,path] = loadData(path)
if nargin < 1
    [fname,pname] = uigetfile('*.*');
    path = [pname,filesep,fname];
end
    
info = imfinfo(path);
w = info(1).Width; h = info(1).Height;
bd = info(1).BitDepth/8;
nf = floor(info(1).FileSize/(bd*w*h));
keepReading = 1; k = 1;
im = zeros(h,w,nf);
warning('off')
while keepReading 
    try
        im(:,:,k) = imread(path,k);
        k = k+1;
    catch
        keepReading = 0;
        disp('Finished reading... ')
        disp(['Stack size : ',num2str(size(im))])
    end
end
warning('on')