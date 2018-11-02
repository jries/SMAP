function out=myxcorr(in1,in2)
out=0*in1;
li=length(in1);
for k=1:length(in1)

    out(k)=sum(in1(1:li-k+1).*in2(k:li));
end