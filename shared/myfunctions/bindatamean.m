function [yy,ss]=bindatamean(x,y,xx,w)

if nargin<4 || isempty(w)
    w=ones(size(x));
end
%ss SEM 

badind=isnan(x)|isnan(y)|isinf(x)|isinf(y);
x(badind)=[];
y(badind)=[];
x(end+1)=inf;
ind1=1;
xn=(xx(1:end-1)+xx(2:end))/2;
xh(1)=-inf; xh(2:length(xn)+1)=xn;xh(end+1)=inf;
xn=xh;
ind2=1;
yy=zeros(size(xx));
ss=zeros(size(xx));
for k=1:length(xx)
    while x(ind2)<xn(k+1) && ind2<length(x)
        ind2=ind2+1;
    end
    indm=ind1:ind2-1;
%     yy(k)=mean(y(indm));
    if ~isempty(indm)
        yy(k)=sum(y(indm).*w(indm))/sum(w(indm));
        ss(k)=sqrt(sum(w(indm).^2.*var(y(indm)))/sum(w(indm))^2);
    else
        ss(k)=0;
    end
    ind1=ind2;
    if ind2>=length(x)
        break
    end

end

    