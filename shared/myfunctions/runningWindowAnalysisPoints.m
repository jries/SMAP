function y=runningWindowAnalysisPoints(x,windowsize,mode)
if nargin<3
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
    case 'robustmean'
        fh=@robustMean;
    otherwise
        disp('bindata.m: mode not known. try parameter as function handle.')
        fh=mode;
end

% [xs,sind]=sort(x);
% x=xs;
% y=y(sind);

badind=isnan(x)|isinf(x);
x(badind)=[];

ws2=floor(windowsize/2);
y=zeros(size(x));
for k=1:length(x)
    ind1=max(1,k-ws2);
    ind2=min(k+ws2,length(x));
    y(k)=fh(x(ind1:ind2));
end

% xn=xx+windowsize/2;
% xn=(xx(1:end-1)+xx(2:end))/2;
% xn(2:end+1)=xn;
% xn(1)=-inf; xn(end+1)=inf;

% for k=1:length(xx)
%     ind=x>=xx(k)-windowsize/2&x<xx(k)+windowsize/2;
%     try
%     yy(k)=fh(y(ind));
%     catch
%         yy(k)=NaN;
%     end
% end
% 