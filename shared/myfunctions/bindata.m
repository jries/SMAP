function [yy,xwin]=bindata(x,y,xx,mode)
if nargin<4
    mode='mean';
end
switch mode
    case 'sum'
        fh=@sum;
    case 'mean'
        fh=@mean;
    case 'median'
        fh=@median;
    case 'geomean'
        fh='geomean';
    case 'std'
        fh=@std;
    case 'robustmean'
        fh=@robustMean;
    case 'rms'
        fh=@rmshere;
    case 'rootsumsquare'
        fh=@rootsumsquarehere;
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
xn(2:end+1)=xn;
xn(1)=-inf; xn(end+1)=inf;
yy=zeros(size(xx));
for k=1:length(xx)
    ind=x>=xn(k)&x<xn(k+1);
    try
    yy(k)=fh(y(ind));
    catch
        yy(k)=NaN;
    end
end
xwin=diff(xn);
end

function out=rmshere(in)
out=sqrt(mean(in.^2));
end

function out=rootsumsquarehere(in)
out=sqrt(sum(in.^2));
end
    