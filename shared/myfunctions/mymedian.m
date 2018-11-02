function out=mymedian(in,dim,maxnum)
shiftin=shiftdim(in,dim-1);
s=size(shiftin);
if length(s)>2
    out=zeros(s(2),s(3),class(in));
for k=1:s(3)
    m=fast_median(shiftin(:,:,k));
    out(:,k)=m;
end
else
    out=in;
end
