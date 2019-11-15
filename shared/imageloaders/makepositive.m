function out=makepositive(in)
if isa(in,'int16')
    out=single(in);
    out(out<0)=out(out<0)+2^16;
else
    out=in;
end
end