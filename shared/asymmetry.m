function [out,asymmr,dir]=asymmetry(imin,subtractmin)
if nargin<2
    subtractmin=false;
end
if isempty(imin)
    out=0;dir=0;asymmr=0;
    return
end
if iscell(imin)
    for k=numel(imin):-1:1
        [out{k},asymmr{k},dir{k}]=asymmetryi(imin{k},subtractmin);
    end
    out=reshape(out,size(imin));dir=reshape(dir,size(imin));
else
    [out,asymmr,dir]=asymmetryi(imin,subtractmin);
end
end
function [out,asymmr,dir]=asymmetryi(imin,subtractmin)
    if (ndims(imin)>2)
        for k=size(imin,3):-1:1
            [out(k),asymmr(k),dir(k)]=asymmetryi(imin(:,:,k),subtractmin);
        end
    else
        [out,asymmr,dir]=asymmetryi2(imin,subtractmin);
    end
end
function [out,asymmr,dir]=asymmetryi2(imin,subtractmin)
persistent X2 Y2 XY sx
s=size(imin);
if isempty(X2) || any(sx~=s)
    n=(0:s(1)-1)-s(1)/2+.5;
    [X,Y]=meshgrid(n,n);
    Xs=(X+Y)/2;
    Ys=(-X+Y)/2;
    X2=Xs.*Xs;
    Y2=Ys.*Ys;
    XY=Xs.*Ys;
    sx=size(X);
end
if subtractmin
    imin=imin-max(0,min(imin(:)));
end
a11=sum(sum(X2.*imin));
a22=sum(sum(Y2.*imin));
a21=sum(sum(XY.*imin));
a=[a11,a21;a21,a22];
[e,v]=eig(a);
if (v(2,2)+v(1,1))>0
    out=abs((v(2,2)-v(1,1))/(v(2,2)+v(1,1)));
else
    out=0;
end
asymmr=abs((a22-a11)/(a22+a11));
    dir=mod(atan2(e(2,2),e(1,2)),pi)/pi*180;
% end
end