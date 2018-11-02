function [out]=detbestlength2(FFTrv,FFTiv,IFFTiv,size1,size2,isreal1,isreal2)
% [out]=detbestlength2(FFTrv,FFTiv,IFFTiv,size1,size2,isreal1,isreal2)
% Determine the best parameters for Overlap-Add FFT-based convolution.
%
% INPUT
% FFTrv:   vector with costs of FFT for real 1d vectors
% FFTiv:   vector with costs of FFT for complex 1d vectors
% IFFTiv:  vector with costs of IFFT for complex 1d vectors
% size1:   size(first_image)
% size2:   size(second_image)
% isreal1: 1 if first image is real, 0 otherwise (complex)
% isreal2: 1 if second image is real, 0 otherwise (complex)
% OUTPUT
% out:    the optimized parameters:
%         out.inverse:     if 1 the two input have to be inverted
%         out.fftxfirst:   if one the image has to be fft first along
%                          x-dimension
%         out.ifftxfirst:  if one the product of spectra has to be ifft
%                          first along x-dimensio
%         out.nfftx:       the best length for fft transform along
%                          x-dimension
%         out.nffty:       the best length for fft transform along
%                          y-dimension
%         out.filterxfirst if 1 the filter has to be fft fisrt alng
%                          x-dimension
%

out           = [];
% the 3 input vectors have to be the same length
L             = length(FFTrv);
% a default value (just as Inf)
infinitevalue = 99*10^99;
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%----------------------------------------------------- a image and b filter
if isreal1 && isreal2
    ax = size1(1);
    ay = size1(2);
    bx = size2(1);
    by = size2(2);

    val0 = infinitevalue;

    for ii=1:L
        for jj=1:L
            if FFTrv(ii)~=0 && FFTrv(jj)~=0
                Lx    = ii-bx+1;
                Ly    = jj-by+1;
                if Lx>0 && Ly>0
                    nx    = ceil(ax/Lx);
                    ny    = ceil(ay/Ly);

                    cv1 = nx*ny*(ii*FFTrv(jj) + jj*FFTiv(ii) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (ii*FFTrv(jj) + jj*FFTiv(ii));
                    cv2 = nx*ny*(ii*FFTrv(jj) + jj*FFTiv(ii) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (ii*FFTrv(jj) + jj*FFTiv(ii));
                    cv3 = nx*ny*(jj*FFTrv(ii) + ii*FFTiv(jj) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (ii*FFTrv(jj) + jj*FFTiv(ii));
                    cv4 = nx*ny*(jj*FFTrv(ii) + ii*FFTiv(jj) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (ii*FFTrv(jj) + jj*FFTiv(ii));

                    cv5 = nx*ny*(ii*FFTrv(jj) + jj*FFTiv(ii) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (jj*FFTrv(ii) + ii*FFTiv(jj));
                    cv6 = nx*ny*(ii*FFTrv(jj) + jj*FFTiv(ii) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (jj*FFTrv(ii) + ii*FFTiv(jj));
                    cv7 = nx*ny*(jj*FFTrv(ii) + ii*FFTiv(jj) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (jj*FFTrv(ii) + ii*FFTiv(jj));
                    cv8 = nx*ny*(jj*FFTrv(ii) + ii*FFTiv(jj) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (jj*FFTrv(ii) + ii*FFTiv(jj));

                    if cv1<val0
                        val0 = cv1;
                        x    = ii;
                        y    = jj;
                        z    = 1;
                        t    = 1;
                    end
                    if cv2<val0
                        val0 = cv2;
                        x    = ii;
                        y    = jj;
                        z    = 2;
                        t    = 1;
                    end
                    if cv3<val0
                        val0 = cv3;
                        x    = ii;
                        y    = jj;
                        z    = 3;
                        t    = 1;
                    end
                    if cv4<val0
                        val0 = cv4;
                        x    = ii;
                        y    = jj;
                        z    = 4;
                        t    = 1;
                    end
                    if cv5<val0
                        val0 = cv5;
                        x    = ii;
                        y    = jj;
                        z    = 1;
                        t    = 0;
                    end
                    if cv6<val0
                        val0 = cv6;
                        x    = ii;
                        y    = jj;
                        z    = 2;
                        t    = 0;
                    end
                    if cv7<val0
                        val0 = cv7;
                        x    = ii;
                        y    = jj;
                        z    = 3;
                        t    = 0;
                    end
                    if cv8<val0
                        val0 = cv8;
                        x    = ii;
                        y    = jj;
                        z    = 4;
                        t    = 0;
                    end
                end
            end
        end
    end
    %----------------------------------------------------- a filter and b image
    ax = size2(1);
    ay = size2(2);
    bx = size1(1);
    by = size1(2);


    val1 = infinitevalue;

    for ii=1:L
        for jj=1:L
            if FFTrv(ii)~=0 && FFTrv(jj)~=0
                Lx    = ii-bx+1;
                Ly    = jj-by+1;
                if Lx>0 && Ly>0
                    nx    = ceil(ax/Lx);
                    ny    = ceil(ay/Ly);

                    cv1 = nx*ny*(ii*FFTrv(jj) + jj*FFTiv(ii) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (ii*FFTrv(jj) + jj*FFTiv(ii));
                    cv2 = nx*ny*(ii*FFTrv(jj) + jj*FFTiv(ii) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (ii*FFTrv(jj) + jj*FFTiv(ii));
                    cv3 = nx*ny*(jj*FFTrv(ii) + ii*FFTiv(jj) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (ii*FFTrv(jj) + jj*FFTiv(ii));
                    cv4 = nx*ny*(jj*FFTrv(ii) + ii*FFTiv(jj) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (ii*FFTrv(jj) + jj*FFTiv(ii));

                    cv5 = nx*ny*(ii*FFTrv(jj) + jj*FFTiv(ii) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (jj*FFTrv(ii) + ii*FFTiv(jj));
                    cv6 = nx*ny*(ii*FFTrv(jj) + jj*FFTiv(ii) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (jj*FFTrv(ii) + ii*FFTiv(jj));
                    cv7 = nx*ny*(jj*FFTrv(ii) + ii*FFTiv(jj) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (jj*FFTrv(ii) + ii*FFTiv(jj));
                    cv8 = nx*ny*(jj*FFTrv(ii) + ii*FFTiv(jj) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (jj*FFTrv(ii) + ii*FFTiv(jj));

                    if cv1<val1
                        val1 = cv1;
                        x    = ii;
                        y    = jj;
                        z    = 1;
                        t    = 1;
                    end
                    if cv2<val1
                        val1 = cv2;
                        x    = ii;
                        y    = jj;
                        z    = 2;
                        t    = 1;
                    end
                    if cv3<val1
                        val1 = cv3;
                        x    = ii;
                        y    = jj;
                        z    = 3;
                        t    = 1;
                    end
                    if cv4<val1
                        val1 = cv4;
                        x    = ii;
                        y    = jj;
                        z    = 4;
                        t    = 1;
                    end
                    if cv5<val1
                        val1 = cv5;
                        x    = ii;
                        y    = jj;
                        z    = 1;
                        t    = 0;
                    end
                    if cv6<val1
                        val1 = cv6;
                        x    = ii;
                        y    = jj;
                        z    = 2;
                        t    = 0;
                    end
                    if cv7<val1
                        val1 = cv7;
                        x    = ii;
                        y    = jj;
                        z    = 3;
                        t    = 0;
                    end
                    if cv8<val1
                        val1 = cv8;
                        x    = ii;
                        y    = jj;
                        z    = 4;
                        t    = 0;
                    end
                end
            end
        end
    end
    %--------------------------------------------------------------------------
    if val1<val0
        out.inverse = 1;
        if z==1 || z==2
            out.fftxfirst = 0;
        else
            out.fftxfirst = 1;
        end
        if z==1 || z==3
            out.ifftxfirst = 1;
        else
            out.ifftxfirst = 0;
        end
    else
        out.inverse = 0;
        if z==1 || z==2
            out.fftxfirst = 0;
        else
            out.fftxfirst = 1;
        end
        if z==1 || z==3
            out.ifftxfirst = 1;
        else
            out.ifftxfirst = 0;
        end
    end
    out.nfftx = x;
    out.nffty = y;
    if t==1
        out.filterxfirst = 0;
    else
        out.filterxfirst = 1;
    end
    return;
end
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%----------------------------------------------------- a image and b filter
if ~isreal1 && ~isreal2
    ax = size1(1);
    ay = size1(2);
    bx = size2(1);
    by = size2(2);

    val0 = infinitevalue;

    for ii=1:L
        for jj=1:L
            if FFTrv(ii)~=0 && FFTrv(jj)~=0
                Lx    = ii-bx+1;
                Ly    = jj-by+1;
                if Lx>0 && Ly>0
                    nx    = ceil(ax/Lx);
                    ny    = ceil(ay/Ly);

                    cv1 = nx*ny*(ii*FFTiv(jj) + jj*FFTiv(ii) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (ii*FFTiv(jj) + jj*FFTiv(ii));
                    cv2 = nx*ny*(ii*FFTiv(jj) + jj*FFTiv(ii) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (ii*FFTiv(jj) + jj*FFTiv(ii));
                    cv3 = nx*ny*(jj*FFTiv(ii) + ii*FFTiv(jj) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (ii*FFTiv(jj) + jj*FFTiv(ii));
                    cv4 = nx*ny*(jj*FFTiv(ii) + ii*FFTiv(jj) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (ii*FFTiv(jj) + jj*FFTiv(ii));

                    cv5 = nx*ny*(ii*FFTiv(jj) + jj*FFTiv(ii) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (jj*FFTiv(ii) + ii*FFTiv(jj));
                    cv6 = nx*ny*(ii*FFTiv(jj) + jj*FFTiv(ii) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (jj*FFTiv(ii) + ii*FFTiv(jj));
                    cv7 = nx*ny*(jj*FFTiv(ii) + ii*FFTiv(jj) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (jj*FFTiv(ii) + ii*FFTiv(jj));
                    cv8 = nx*ny*(jj*FFTiv(ii) + ii*FFTiv(jj) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (jj*FFTiv(ii) + ii*FFTiv(jj));

                    if cv1<val0
                        val0 = cv1;
                        x    = ii;
                        y    = jj;
                        z    = 1;
                        t    = 1;
                    end
                    if cv2<val0
                        val0 = cv2;
                        x    = ii;
                        y    = jj;
                        z    = 2;
                        t    = 1;
                    end
                    if cv3<val0
                        val0 = cv3;
                        x    = ii;
                        y    = jj;
                        z    = 3;
                        t    = 1;
                    end
                    if cv4<val0
                        val0 = cv4;
                        x    = ii;
                        y    = jj;
                        z    = 4;
                        t    = 1;
                    end
                    if cv5<val0
                        val0 = cv5;
                        x    = ii;
                        y    = jj;
                        z    = 1;
                        t    = 0;
                    end
                    if cv6<val0
                        val0 = cv6;
                        x    = ii;
                        y    = jj;
                        z    = 2;
                        t    = 0;
                    end
                    if cv7<val0
                        val0 = cv7;
                        x    = ii;
                        y    = jj;
                        z    = 3;
                        t    = 0;
                    end
                    if cv8<val0
                        val0 = cv8;
                        x    = ii;
                        y    = jj;
                        z    = 4;
                        t    = 0;
                    end
                end
            end
        end
    end
    %----------------------------------------------------- a filter and b image
    ax = size2(1);
    ay = size2(2);
    bx = size1(1);
    by = size1(2);


    val1 = infinitevalue;

    for ii=1:L
        for jj=1:L
            if FFTiv(ii)~=0 && FFTiv(jj)~=0
                Lx    = ii-bx+1;
                Ly    = jj-by+1;
                if Lx>0 && Ly>0
                    nx    = ceil(ax/Lx);
                    ny    = ceil(ay/Ly);

                    cv1 = nx*ny*(ii*FFTiv(jj) + jj*FFTiv(ii) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (ii*FFTiv(jj) + jj*FFTiv(ii));
                    cv2 = nx*ny*(ii*FFTiv(jj) + jj*FFTiv(ii) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (ii*FFTiv(jj) + jj*FFTiv(ii));
                    cv3 = nx*ny*(jj*FFTiv(ii) + ii*FFTiv(jj) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (ii*FFTiv(jj) + jj*FFTiv(ii));
                    cv4 = nx*ny*(jj*FFTiv(ii) + ii*FFTiv(jj) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (ii*FFTiv(jj) + jj*FFTiv(ii));

                    cv5 = nx*ny*(ii*FFTiv(jj) + jj*FFTiv(ii) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (jj*FFTiv(ii) + ii*FFTiv(jj));
                    cv6 = nx*ny*(ii*FFTiv(jj) + jj*FFTiv(ii) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (jj*FFTiv(ii) + ii*FFTiv(jj));
                    cv7 = nx*ny*(jj*FFTiv(ii) + ii*FFTiv(jj) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (jj*FFTiv(ii) + ii*FFTiv(jj));
                    cv8 = nx*ny*(jj*FFTiv(ii) + ii*FFTiv(jj) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (jj*FFTiv(ii) + ii*FFTiv(jj));

                    if cv1<val1
                        val1 = cv1;
                        x    = ii;
                        y    = jj;
                        z    = 1;
                        t    = 1;
                    end
                    if cv2<val1
                        val1 = cv2;
                        x    = ii;
                        y    = jj;
                        z    = 2;
                        t    = 1;
                    end
                    if cv3<val1
                        val1 = cv3;
                        x    = ii;
                        y    = jj;
                        z    = 3;
                        t    = 1;
                    end
                    if cv4<val1
                        val1 = cv4;
                        x    = ii;
                        y    = jj;
                        z    = 4;
                        t    = 1;
                    end
                    if cv5<val1
                        val1 = cv5;
                        x    = ii;
                        y    = jj;
                        z    = 1;
                        t    = 0;
                    end
                    if cv6<val1
                        val1 = cv6;
                        x    = ii;
                        y    = jj;
                        z    = 2;
                        t    = 0;
                    end
                    if cv7<val1
                        val1 = cv7;
                        x    = ii;
                        y    = jj;
                        z    = 3;
                        t    = 0;
                    end
                    if cv8<val1
                        val1 = cv8;
                        x    = ii;
                        y    = jj;
                        z    = 4;
                        t    = 0;
                    end
                end
            end
        end
    end
    %--------------------------------------------------------------------------
    if val1<val0
        out.inverse = 1;
        if z==1 || z==2
            out.fftxfirst = 0;
        else
            out.fftxfirst = 1;
        end
        if z==1 || z==3
            out.ifftxfirst = 1;
        else
            out.ifftxfirst = 0;
        end
    else
        out.inverse = 0;
        if z==1 || z==2
            out.fftxfirst = 0;
        else
            out.fftxfirst = 1;
        end
        if z==1 || z==3
            out.ifftxfirst = 1;
        else
            out.ifftxfirst = 0;
        end
    end
    out.nfftx = x;
    out.nffty = y;
    if t==1
        out.filterxfirst = 0;
    else
        out.filterxfirst = 1;
    end
    return;
end
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%----------------------------------------------------- a image and b filter
if isreal1 && ~isreal2
    ax = size1(1);
    ay = size1(2);
    bx = size2(1);
    by = size2(2);

    val0 = infinitevalue;

    for ii=1:L
        for jj=1:L
            if FFTrv(ii)~=0 && FFTrv(jj)~=0
                Lx    = ii-bx+1;
                Ly    = jj-by+1;
                if Lx>0 && Ly>0
                    nx    = ceil(ax/Lx);
                    ny    = ceil(ay/Ly);

                    cv1 = nx*ny*(ii*FFTrv(jj) + jj*FFTiv(ii) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (ii*FFTiv(jj) + jj*FFTiv(ii));
                    cv2 = nx*ny*(ii*FFTrv(jj) + jj*FFTiv(ii) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (ii*FFTiv(jj) + jj*FFTiv(ii));
                    cv3 = nx*ny*(jj*FFTrv(ii) + ii*FFTiv(jj) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (ii*FFTiv(jj) + jj*FFTiv(ii));
                    cv4 = nx*ny*(jj*FFTrv(ii) + ii*FFTiv(jj) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (ii*FFTiv(jj) + jj*FFTiv(ii));

                    cv5 = nx*ny*(ii*FFTrv(jj) + jj*FFTiv(ii) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (jj*FFTiv(ii) + ii*FFTiv(jj));
                    cv6 = nx*ny*(ii*FFTrv(jj) + jj*FFTiv(ii) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (jj*FFTiv(ii) + ii*FFTiv(jj));
                    cv7 = nx*ny*(jj*FFTrv(ii) + ii*FFTiv(jj) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (jj*FFTiv(ii) + ii*FFTiv(jj));
                    cv8 = nx*ny*(jj*FFTrv(ii) + ii*FFTiv(jj) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (jj*FFTiv(ii) + ii*FFTiv(jj));

                    if cv1<val0
                        val0 = cv1;
                        x    = ii;
                        y    = jj;
                        z    = 1;
                        t    = 1;
                    end
                    if cv2<val0
                        val0 = cv2;
                        x    = ii;
                        y    = jj;
                        z    = 2;
                        t    = 1;
                    end
                    if cv3<val0
                        val0 = cv3;
                        x    = ii;
                        y    = jj;
                        z    = 3;
                        t    = 1;
                    end
                    if cv4<val0
                        val0 = cv4;
                        x    = ii;
                        y    = jj;
                        z    = 4;
                        t    = 1;
                    end
                    if cv5<val0
                        val0 = cv5;
                        x    = ii;
                        y    = jj;
                        z    = 1;
                        t    = 0;
                    end
                    if cv6<val0
                        val0 = cv6;
                        x    = ii;
                        y    = jj;
                        z    = 2;
                        t    = 0;
                    end
                    if cv7<val0
                        val0 = cv7;
                        x    = ii;
                        y    = jj;
                        z    = 3;
                        t    = 0;
                    end
                    if cv8<val0
                        val0 = cv8;
                        x    = ii;
                        y    = jj;
                        z    = 4;
                        t    = 0;
                    end
                end
            end
        end
    end
    %----------------------------------------------------- a filter and b image
    ax = size2(1);
    ay = size2(2);
    bx = size1(1);
    by = size1(2);


    val1 = infinitevalue;

    for ii=1:L
        for jj=1:L
            if FFTrv(ii)~=0 && FFTrv(jj)~=0
                Lx    = ii-bx+1;
                Ly    = jj-by+1;
                if Lx>0 && Ly>0
                    nx    = ceil(ax/Lx);
                    ny    = ceil(ay/Ly);

                    cv1 = nx*ny*(ii*FFTiv(jj) + jj*FFTiv(ii) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (ii*FFTrv(jj) + jj*FFTiv(ii));
                    cv2 = nx*ny*(ii*FFTiv(jj) + jj*FFTiv(ii) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (ii*FFTrv(jj) + jj*FFTiv(ii));
                    cv3 = nx*ny*(jj*FFTiv(ii) + ii*FFTiv(jj) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (ii*FFTrv(jj) + jj*FFTiv(ii));
                    cv4 = nx*ny*(jj*FFTiv(ii) + ii*FFTiv(jj) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (ii*FFTrv(jj) + jj*FFTiv(ii));

                    cv5 = nx*ny*(ii*FFTiv(jj) + jj*FFTiv(ii) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (jj*FFTrv(ii) + ii*FFTiv(jj));
                    cv6 = nx*ny*(ii*FFTiv(jj) + jj*FFTiv(ii) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (jj*FFTrv(ii) + ii*FFTiv(jj));
                    cv7 = nx*ny*(jj*FFTiv(ii) + ii*FFTiv(jj) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (jj*FFTrv(ii) + ii*FFTiv(jj));
                    cv8 = nx*ny*(jj*FFTiv(ii) + ii*FFTiv(jj) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (jj*FFTrv(ii) + ii*FFTiv(jj));

                    if cv1<val1
                        val1 = cv1;
                        x    = ii;
                        y    = jj;
                        z    = 1;
                        t    = 1;
                    end
                    if cv2<val1
                        val1 = cv2;
                        x    = ii;
                        y    = jj;
                        z    = 2;
                        t    = 1;
                    end
                    if cv3<val1
                        val1 = cv3;
                        x    = ii;
                        y    = jj;
                        z    = 3;
                        t    = 1;
                    end
                    if cv4<val1
                        val1 = cv4;
                        x    = ii;
                        y    = jj;
                        z    = 4;
                        t    = 1;
                    end
                    if cv5<val1
                        val1 = cv5;
                        x    = ii;
                        y    = jj;
                        z    = 1;
                        t    = 0;
                    end
                    if cv6<val1
                        val1 = cv6;
                        x    = ii;
                        y    = jj;
                        z    = 2;
                        t    = 0;
                    end
                    if cv7<val1
                        val1 = cv7;
                        x    = ii;
                        y    = jj;
                        z    = 3;
                        t    = 0;
                    end
                    if cv8<val1
                        val1 = cv8;
                        x    = ii;
                        y    = jj;
                        z    = 4;
                        t    = 0;
                    end
                end
            end
        end
    end
    %--------------------------------------------------------------------------
    if val1<val0
        out.inverse = 1;
        if z==1 || z==2
            out.fftxfirst = 0;
        else
            out.fftxfirst = 1;
        end
        if z==1 || z==3
            out.ifftxfirst = 1;
        else
            out.ifftxfirst = 0;
        end
    else
        out.inverse = 0;
        if z==1 || z==2
            out.fftxfirst = 0;
        else
            out.fftxfirst = 1;
        end
        if z==1 || z==3
            out.ifftxfirst = 1;
        else
            out.ifftxfirst = 0;
        end
    end
    out.nfftx = x;
    out.nffty = y;
    if t==1
        out.filterxfirst = 0;
    else
        out.filterxfirst = 1;
    end
    return;
end
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%----------------------------------------------------- a image and b filter
if ~isreal1 && isreal2
    ax = size1(1);
    ay = size1(2);
    bx = size2(1);
    by = size2(2);

    val0 = infinitevalue;

    for ii=1:L
        for jj=1:L
            if FFTrv(ii)~=0 && FFTrv(jj)~=0
                Lx    = ii-bx+1;
                Ly    = jj-by+1;
                if Lx>0 && Ly>0
                    nx    = ceil(ax/Lx);
                    ny    = ceil(ay/Ly);

                    cv1 = nx*ny*(ii*FFTiv(jj) + jj*FFTiv(ii) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (ii*FFTrv(jj) + jj*FFTiv(ii));
                    cv2 = nx*ny*(ii*FFTiv(jj) + jj*FFTiv(ii) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (ii*FFTrv(jj) + jj*FFTiv(ii));
                    cv3 = nx*ny*(jj*FFTiv(ii) + ii*FFTiv(jj) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (ii*FFTrv(jj) + jj*FFTiv(ii));
                    cv4 = nx*ny*(jj*FFTiv(ii) + ii*FFTiv(jj) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (ii*FFTrv(jj) + jj*FFTiv(ii));

                    cv5 = nx*ny*(ii*FFTiv(jj) + jj*FFTiv(ii) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (jj*FFTrv(ii) + ii*FFTiv(jj));
                    cv6 = nx*ny*(ii*FFTiv(jj) + jj*FFTiv(ii) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (jj*FFTrv(ii) + ii*FFTiv(jj));
                    cv7 = nx*ny*(jj*FFTiv(ii) + ii*FFTiv(jj) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (jj*FFTrv(ii) + ii*FFTiv(jj));
                    cv8 = nx*ny*(jj*FFTiv(ii) + ii*FFTiv(jj) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (jj*FFTrv(ii) + ii*FFTiv(jj));

                    if cv1<val0
                        val0 = cv1;
                        x    = ii;
                        y    = jj;
                        z    = 1;
                        t    = 1;
                    end
                    if cv2<val0
                        val0 = cv2;
                        x    = ii;
                        y    = jj;
                        z    = 2;
                        t    = 1;
                    end
                    if cv3<val0
                        val0 = cv3;
                        x    = ii;
                        y    = jj;
                        z    = 3;
                        t    = 1;
                    end
                    if cv4<val0
                        val0 = cv4;
                        x    = ii;
                        y    = jj;
                        z    = 4;
                        t    = 1;
                    end
                    if cv5<val0
                        val0 = cv5;
                        x    = ii;
                        y    = jj;
                        z    = 1;
                        t    = 0;
                    end
                    if cv6<val0
                        val0 = cv6;
                        x    = ii;
                        y    = jj;
                        z    = 2;
                        t    = 0;
                    end
                    if cv7<val0
                        val0 = cv7;
                        x    = ii;
                        y    = jj;
                        z    = 3;
                        t    = 0;
                    end
                    if cv8<val0
                        val0 = cv8;
                        x    = ii;
                        y    = jj;
                        z    = 4;
                        t    = 0;
                    end
                end
            end
        end
    end
    %----------------------------------------------------- a filter and b image
    ax = size2(1);
    ay = size2(2);
    bx = size1(1);
    by = size1(2);


    val1 = infinitevalue;

    for ii=1:L
        for jj=1:L
            if FFTrv(ii)~=0 && FFTrv(jj)~=0
                Lx    = ii-bx+1;
                Ly    = jj-by+1;
                if Lx>0 && Ly>0
                    nx    = ceil(ax/Lx);
                    ny    = ceil(ay/Ly);

                    cv1 = nx*ny*(ii*FFTrv(jj) + jj*FFTiv(ii) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (ii*FFTiv(jj) + jj*FFTiv(ii));
                    cv2 = nx*ny*(ii*FFTrv(jj) + jj*FFTiv(ii) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (ii*FFTiv(jj) + jj*FFTiv(ii));
                    cv3 = nx*ny*(jj*FFTrv(ii) + ii*FFTiv(jj) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (ii*FFTiv(jj) + jj*FFTiv(ii));
                    cv4 = nx*ny*(jj*FFTrv(ii) + ii*FFTiv(jj) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (ii*FFTiv(jj) + jj*FFTiv(ii));

                    cv5 = nx*ny*(ii*FFTrv(jj) + jj*FFTiv(ii) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (jj*FFTiv(ii) + ii*FFTiv(jj));
                    cv6 = nx*ny*(ii*FFTrv(jj) + jj*FFTiv(ii) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (jj*FFTiv(ii) + ii*FFTiv(jj));
                    cv7 = nx*ny*(jj*FFTrv(ii) + ii*FFTiv(jj) + jj*IFFTiv(ii) + ii*IFFTiv(jj)) + (jj*FFTiv(ii) + ii*FFTiv(jj));
                    cv8 = nx*ny*(jj*FFTrv(ii) + ii*FFTiv(jj) + ii*IFFTiv(jj) + jj*IFFTiv(ii)) + (jj*FFTiv(ii) + ii*FFTiv(jj));

                    if cv1<val1
                        val1 = cv1;
                        x    = ii;
                        y    = jj;
                        z    = 1;
                        t    = 1;
                    end
                    if cv2<val1
                        val1 = cv2;
                        x    = ii;
                        y    = jj;
                        z    = 2;
                        t    = 1;
                    end
                    if cv3<val1
                        val1 = cv3;
                        x    = ii;
                        y    = jj;
                        z    = 3;
                        t    = 1;
                    end
                    if cv4<val1
                        val1 = cv4;
                        x    = ii;
                        y    = jj;
                        z    = 4;
                        t    = 1;
                    end
                    if cv5<val1
                        val1 = cv5;
                        x    = ii;
                        y    = jj;
                        z    = 1;
                        t    = 0;
                    end
                    if cv6<val1
                        val1 = cv6;
                        x    = ii;
                        y    = jj;
                        z    = 2;
                        t    = 0;
                    end
                    if cv7<val1
                        val1 = cv7;
                        x    = ii;
                        y    = jj;
                        z    = 3;
                        t    = 0;
                    end
                    if cv8<val1
                        val1 = cv8;
                        x    = ii;
                        y    = jj;
                        z    = 4;
                        t    = 0;
                    end
                end
            end
        end
    end
    %--------------------------------------------------------------------------
    if val1<val0
        out.inverse = 1;
        if z==1 || z==2
            out.fftxfirst = 0;
        else
            out.fftxfirst = 1;
        end
        if z==1 || z==3
            out.ifftxfirst = 1;
        else
            out.ifftxfirst = 0;
        end
    else
        out.inverse = 0;
        if z==1 || z==2
            out.fftxfirst = 0;
        else
            out.fftxfirst = 1;
        end
        if z==1 || z==3
            out.ifftxfirst = 1;
        else
            out.ifftxfirst = 0;
        end
    end
    out.nfftx = x;
    out.nffty = y;
    if t==1
        out.filterxfirst = 0;
    else
        out.filterxfirst = 1;
    end
    return;
end
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
