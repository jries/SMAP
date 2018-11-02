function n=myquantilefast(in,q,maxnum)
if isempty(in)
    n=[];
    return
end
if nargin>2 && maxnum<numel(in)*.75
    dn=ceil(numel(in)/maxnum);
    in2=in(1:dn:end);
else
in2=in;
end
nin=numel(in2);
n=zeros(size(q));
for k=1:length(q)
    p=(ceil(nin*q(k)));
    d=nth_element(in2(:),p);
    n(k)=d(p);
end
