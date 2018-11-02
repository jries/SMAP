function [out, meanv, stdv]=bintrace(in,win)
l=floor(length(in)/win);
out=zeros(l,1);meanv=zeros(l,1);stdv=zeros(l,1);
for i=1:l
    out(i)=sum(in((i-1)*win+1:i*win));
    meanv(i)=mean(in((i-1)*win+1:i*win));
    stdv(i)=std(in((i-1)*win+1:i*win));
end