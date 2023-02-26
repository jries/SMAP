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

if ~issorted(x)
    [xs,sind]=sort(x);
    x=xs;
    y=y(sind);
end

badind=isnan(x)|isnan(y)|isinf(x)|isinf(y);
x(badind)=[];
y(badind)=[];

xn=(xx(1:end-1)+xx(2:end))/2;
xn(2:end+1)=xn;
xn(1)=-inf; xn(end+1)=inf;
yy=zeros(size(xx));

if length(x)<1e5 %old way

    for k=1:length(xx)
        ind=x>=xn(k)&x<xn(k+1);
        try
        yy(k)=fh(y(ind));
        catch
            yy(k)=NaN;
        end
    end
else
    ind1=1;
    for k=1:length(xx)
        while ind1<=length(x) && x(ind1)<xn(k)
            ind1=ind1+1;
        end
        ind2=ind1;
        while ind2<=length(x) && x(ind2)<xn(k+1)
            ind2=ind2+1;
        end
        try
            yy(k)=fh(y(ind1:ind2-1));
        catch
            yy(k)=NaN;
        end
        ind1=ind2;
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
    