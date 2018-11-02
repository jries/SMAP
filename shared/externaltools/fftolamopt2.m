function [out]=fftolamopt2(a,b,opt,shape)
% [out]=fftolamopt2(a,b,siz1,siz2,shape)
%
% Overlap-add method FFT-based 2D convolution
% Example:
%   load fftexecutiontimes;                                                        % load FFTrv, FFTiv and IFFTiv in workspace
%   a   = rand(500,500);                                                           % first image
%   b   = rand(340,220);                                                           % second image
%   opt = detbestlength2(FFTrv,FFTiv,IFFTiv,size(a),size(b),isreal(a),isreal(b));  % optimized parameters
%   y0  = fftolamopt2(a,b,opt);                                                    % equivalent to y0 = conv2(a,b);
%
% INPUT
% a:     first image (2D double matrix)
% b:     second image (2D double matrix)
% opt:   the optimized parameters calculated by detbestlength.m function
%        opt = detbestlength(FFTrv,FFTiv,IFFTiv,size(a),size(b));
% shape: returns a subsection of the 2D convolution with size specified by
%        'shape':
%          'full'  - (default) returns the full 2-D convolution,
%          'same'  - returns the central part of the convolution
%                    that is the same size as A.
%          'valid' - returns only those parts of the convolution
%                    that are computed without the zero-padded
%                    edges. size(C) = [ma-mb+1,na-nb+1] when
%                    all(size(A) >= size(B)), otherwise C is empty.
% See also conv2.
% OUTPUT
% out:   2D convolution of a and b matrices: out = conv2(a,b);


% Original size
[z1x,z1y] = size(a);
[z2x,z2y] = size(b);

% Reverse a and b if necessary
if opt.inverse
    atemp = a;
    a     = b;
    b     = atemp;
end

fftorder  = zeros(2,1);
ifftorder = zeros(2,1);
fftsize   = zeros(2,1);
filterord = zeros(2,1);
filtersiz = zeros(2,1);

if (opt.fftxfirst == 1)
    fftorder(1)  = 1;
    fftorder(2)  = 2;
    fftsize(1)   = opt.nfftx;
    fftsize(2)   = opt.nffty;
else
    fftorder(1)  = 2;
    fftorder(2)  = 1;
    fftsize(1)   = opt.nffty;
    fftsize(2)   = opt.nfftx;
end


if (opt.ifftxfirst == 1)
    ifftorder(1) = 1;
    ifftorder(2) = 2;
else
    ifftorder(1) = 2;
    ifftorder(2) = 1;
end

if opt.filterxfirst==1
    filterord(1) = 1;
    filterord(2) = 2;

    filtersiz(1) = opt.nfftx;
    filtersiz(2) = opt.nffty;
else
    filterord(1) = 2;
    filterord(2) = 1;

    filtersiz(1) = opt.nffty;
    filtersiz(2) = opt.nfftx;
end

siz1          = opt.nfftx;
siz2          = opt.nffty;

[ax,ay]       = size(a);
[bx,by]       = size(b);
dimx          = ax+bx-1;
dimy          = ay+by-1;
nfftx         = siz1;
nffty         = siz2;
Lx            = nfftx-bx+1;
Ly            = nffty-by+1;
B             = fft(fft(b,filtersiz(1),filterord(1)),filtersiz(2),filterord(2));
out           = zeros(dimx,dimy);
x0 = 1;
while x0 <= ax
    x1   = min(x0+Lx-1,ax);
    y0   = 1;
    endx = min(dimx,x0+nfftx-1);
    while y0 <= ay
        y1                   = min(y0+Ly-1,ay);
        endy                 = min(dimy,y0+nffty-1);
        X                    = fft(fft(a(x0:x1,y0:y1),fftsize(1),fftorder(1)),fftsize(2),fftorder(2));
        Y                    = ifft(ifft(X.*B,[],ifftorder(1)),[],ifftorder(2));
        out(x0:endx,y0:endy) = out(x0:endx,y0:endy)+Y(1:(endx-x0+1),1:(endy-y0+1));
        y0                   = y0+Ly;
    end
    x0 = x0+Lx;
end
if isreal(a) && isreal(b)
    out=real(out);
end
if nargin<4 || strcmp(shape,'full')
    return;
end
if strcmp(shape,'valid')
    if ((z1x<z2x)||(z1y<z2y))
        out = [];
    else
        px  = z2x;
        py  = z2y;
        out = out(px:px+z1x-z2x,py:py+z1y-z2y);
    end
    return;
end
if strcmp(shape,'same')
    px  = ((z2x-1)+mod((z2x-1),2))/2;
    py  = ((z2y-1)+mod((z2y-1),2))/2;
    out = out(px+1:px+z1x,py+1:py+z1y);
    return;
end