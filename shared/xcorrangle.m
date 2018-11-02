function out=xcorrangle(in1,in2)
if nargin==1
    in2=in1;
end
out=0*in1;
li=length(in1);
for k=1:length(in1)

    out(k)=sum(in1(1:li-k+1).*in2(k:li))+ sum(in1(li-k+2:li).*in2(1:k-1));
end