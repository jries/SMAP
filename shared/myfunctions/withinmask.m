function ind=withinmask(mask,x,y,z)
s=size(mask);
x=round(x);x(x<1)=1;x(x>s(1))=s(1);
y=round(y);y(y<1)=1;y(y>s(2))=s(2);
if nargin==3
ind=mask(sub2ind(size(mask),round(x),round(y)));
else
    z=round(z);z(z<1)=1;z(z>s(3))=s(3);
    ind=mask(sub2ind(size(mask),round(x),round(y),round(z)));
end