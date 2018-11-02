function yy=bindata(x,y,xx,mode)
if nargin<4
    mode='mean';
end
switch mode
    case 'mean'
        fh=@mean;
    case 'median'
        fh=@median;
    case 'geomean'
        fh='geomean';
    case 'std'
        fh=@std;
    otherwise
        disp('bindata.m: mode not known. try parameter as function handle.')
        fh=mode;
end

% [xs,sind]=sort(x);
% x=xs;
% y=y(sind);

badind=isnan(x)|isnan(y)|isinf(x)|isinf(y);
x(badind)=[];
y(badind)=[];

xn=(xx(1:end-1)+xx(2:end))/2;
xn=[-inf xn inf];
yy=zeros(size(xx));
for k=1:length(xx)
    ind=x>xn(k)&x<xn(k+1);
    yy(k)=fh(y(ind));
end
    