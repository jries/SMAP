function out=xcorr2fft(in1,in2)
if nargin==1
    srec=size(in1);
    nfftexp=2^ceil(log2(max(srec)))*2;
    fi1=fft2(in1,nfftexp,nfftexp); 
    cc=fi1.*conj(fi1);
    out=fftshift(ifft2(cc));
else
    srec=size(in1);
    nfftexp=2^ceil(log2(max(srec)))*2;
    fi1=fft2(in1,nfftexp,nfftexp);
    fi2=fft2(in2,nfftexp,nfftexp);    
    cc=fi1.*conj(fi2);
    out=fftshift(ifft2(cc));
end
