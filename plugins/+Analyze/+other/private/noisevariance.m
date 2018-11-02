% NOISEVARIANEC = Calculates the noise variance and mean intensity for an image
%
% SYNOPSIS:
%   [noisevariance, sumintensity] = noisevariance(in, kthres, AvoidCross)
%
%   kthres: frequencies that should be cut away
%        in the high pass filtering [0, sqrt(dimension of in)]
%        the 'corners' of the rectangular image have freq > 1
%
% DEFAULTS:
%   kthres = .9
%   AvoidCross =1
%
% SEE ALSO
%   pcfo
%
% LITERATURE:
%   R. Heintzmann, P. Relich, R. Nieuwenhuizen, K. Lidke, B. Rieger, 
%   Calibrating photon counts from a single image, submitted
%
% Requires matlab toolbox DIPimage to run (free academic download from www.diplib.org)

% (C) Copyright 2004-2016               
%     All rights reserved               Faculty of Applied Sciences
%                                       Delft University of Technology
%                                       Delft, The Netherlands
%                                       &
%                                       Institute of Physical Chemistry, 
%                                       Friedrich Schiller University,
%                                       Jena, Germany
% Bernd Rieger & Rainer Heintzmann

function [TotalVar, ImgSum] = noisevariance(in, kthres, AvoidCross)
borderfraction=0.05;  % real space fraction of pixels to be used. This can avoid some residual artefacts near the border.
CrossWidth=3;  % Width of pixels +/- to avoid in the cross

if ~isa(in,'dip_image');
    in = mat2im(in);
end
if ndims(in)>2
    in = squeeze(in);
    if ndims(in)>2; error('Input image must be 2D.');end
end
if nargin < 3
    AvoidCross = 1;
end
if nargin <2
    kthres =.9;
else
    if kthres <0 | kthres >sqrt(ndims(in))
        error('kcut in [0,sqrt(dimensions of in)]');
    end
end

inOriginal = in;
in = symmetrize(in); %flip and mirror the input to avoid boundary effects in the DFT
f = ft(in);
m = rr(in)>(kthres*size(in,1)/2);  % mask for high-frequency region
if (AvoidCross)  % introduce the blocked regions anyway
    m(abs(xx(m)) <=CrossWidth)=0;
    m(abs(yy(m)) <=CrossWidth)=0;
end
% make the masks smooth to avoid too extensive spreading in real space. Seems to yield a small advantage for the mean value.
m=real(gaussf(m,[CrossWidth CrossWidth]/4));

fra = sum(m)/prod(size(m));   % fraction of all pixels in high-pass filter

nf = m.*f;   % filtered result
nf=ift(nf);  % Not needed any longer, as the power spectral energy can just as well be measured in Fourier space.

amask=newim(inOriginal);
myborder=round(borderfraction/2*size(amask));
amask(myborder(1):end-myborder(1),myborder(2):end-myborder(2))=1;
amaskBig=symmetrize(amask);

TotalVar = mean(abssqr(nf),amaskBig)/fra;   
ImgSum = mean(inOriginal,amask);

return;


