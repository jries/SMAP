% [kcMap,A0Map] = getLocalDcorr(im,tileSize,tileOverlap,r,Ng,figID)
% ---------------------------------------
%
% Estimate the image local cut-off frequency based on decorrelation analysis 
%
% Inputs:
%  im        	2D image to be analyzed
%  tileSize		Size of a tile in pixels
%  tileOverlap  Spatial overlap between two consecutive tiles
%  r           	Fourier space sampling of the analysis (default: r = linspace(0,1,50)
%  Ng			Number of high-pass filtering (default: Ng = 10)
%  figID		If figID > 1, curves will be plotted in figure(figID)
%
% Outputs:
%  kcMax        Estimated cut-off frequency of the image in normalized frequency for each sectors
%  A0			Amplitude of the local maxima of d0 for each sectors
%  kcGM			Estimated cut-off frequency using Geometric-Mean metric for each sectors
%  d0 			Decorrelation function before high-pass filtering for each sectors
%  d			All decorrelation functions for each sectors
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

function [kcMap,A0Map] = getLocalDcorr(im,tileSize,tileOverlap,r,Ng,figID)

if nargin < 6; figID = 0; end
if nargin < 5; Ng = 10; end
if nargin < 4; r = linspace(0,1,50);end

px = 1:(tileSize-tileOverlap+1):(size(im,2)-tileSize);
py = 1:(tileSize-tileOverlap+1):(size(im,1)-tileSize);
kcMap = zeros(length(py)-1,length(px)-1);
A0Map = kcMap;
for xx = 1:length(px)
	for yy = 1:length(py)
        subIm = im(py(yy):py(yy)+tileSize,px(xx):px(xx)+tileSize,1);
        subIm = subIm(1:size(subIm,1)-not(mod(size(subIm,1),2)),1:size(subIm,2)-not(mod(size(subIm,2),2)));
        [kc,A0] = getDcorr(apodImRect(subIm,20),r,Ng,(figID+1)*(figID>0));
        kcMap(yy,xx) = kc;
        A0Map(yy,xx) = A0;
	end
end

if figID
    figure(figID)
    subplot(121)
        imagesc(kcMap); colorbar; title('kcMap')
    subplot(122)
        imagesc(A0Map); colorbar; title('A0Map')
end
    
    