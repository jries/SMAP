%  Copyright (c)2017 Ries Lab, European Molecular Biology Laboratory,
%  Heidelberg.
%  
%  This file is part of GPUmleFit_LM Fitter.
%  
%  GPUmleFit_LM Fitter is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%  
%  GPUmleFit_LM Fitter is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%  
%  You should have received a copy of the GNU General Public License
%  along with GPUmleFit_LM Fitter.  If not, see <http://www.gnu.org/licenses/>.
%  
%  
%  Additional permission under GNU GPL version 3 section 7

%%
function [out] = simSplinePSF(Npixels,coeff,I,bg,cor)
t=tic;
if (nargin <5)
   error('Minimal usage: simSplinePSF(Npixels,coeff,I,bg,cor)');
end

% if (bg == 0)
%     bg = 10^-10;
% end    

Nfits = size(cor,1);
spline_xsize = size(coeff,1);
spline_ysize = size(coeff,2);
spline_zsize = size(coeff,3);
off = floor(((spline_xsize+1)-Npixels)/2);
data = zeros(Npixels,Npixels,Nfits,'single');


for kk = 1:Nfits
    xcenter = cor(kk,1);
    ycenter = cor(kk,2);
    zcenter = cor(kk,3);
    
    xc = -1*(xcenter - Npixels/2+0.5);
    yc = -1*(ycenter - Npixels/2+0.5);
    zc = zcenter - floor(zcenter);
    
    xstart = floor(xc);
    xc = xc - xstart;
    
    ystart = floor(yc);
    yc = yc - ystart;
    

    zstart = floor(zcenter);
    
   
    [delta_f,delta_dxf,delta_ddxf,delta_dyf,delta_ddyf,delta_dzf,delta_ddzf]=computeDelta3Dj_v2(single(xc),single(yc),single(zc));
    
    for ii = 0:Npixels-1
        for jj = 0:Npixels-1
             temp = fAt3Dj_v2(ii+xstart+off,jj+ystart+off,zstart,spline_xsize,spline_ysize,spline_zsize,delta_f,coeff);
             model = temp*I+bg;
             data(ii+1,jj+1,kk)=model;
        end
    end
    if toc(t)>1
        disp(kk/Nfits)
        t=tic;
    end
    
end
out = (poissrnd(data,Npixels,Npixels,Nfits)); 

