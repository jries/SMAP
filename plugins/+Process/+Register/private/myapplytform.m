function [xo,yo]=myapplytform(x,y,transform)
x=double(x);y=double(y);

%if tform then trnsform directly. Else use own structure

switch transform.refpart
    case 1 %all
        ind1=true(size(x));
    case 2 %left
        ind1=x<=transform.midp;
    case 3%right
        ind1=x>=transform.midp;
    case 4 %top
        ind1=y<=transform.midp;
    case 5% bottom
        ind1=y>=transform.midp;
    otherwise %all
        dips('should not happen')
end
xo=zeros(size(x));
yo=zeros(size(y));
[xo(ind1),yo(ind1)]=transformPointsInverse(transform.T,x(ind1),y(ind1));
% [xo(~ind1),yo(~ind1)]=transformPointsForward(transform.T,x(~ind1),y(~ind1));