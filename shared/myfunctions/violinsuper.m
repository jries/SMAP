function hp=violinsuper(varargin)
% 1. edges values for histogram bins
% 2. array with histograms or cell with raw data
% 3. optional: position of violin plot on x-axis
% Name - value pairs:
%'Parent' axis handle
%'lut' lookup table: 'lines' (default), 'jet', 'hsv' 'parula', or any of the many MATLAB LUTs.;
% 'FaceAlpha': opacity of histogrms
%'LineWidth' : of histogrms
% 'statistics','mean' or 'median'
% 'quantile': length of error bars when median is chosen (quantile). For
% mean the standard deviation is used
%'width': width of the violin


p=parseinput(varargin);
if isempty(p.Parent)
    ax=gca;
else
    ax=p.Parent;
end

if iscell(p.histograms)
    inputhist=false;
    [histograms, edges]=makehistogramsi(p.histograms,p.binpositions);
    dat=p.histograms;
else
    inputhist=true;
    histograms=p.histograms;
    edges=p.binpositions;
end

numhist=size(histograms,2);
h0=[zeros(size(histograms,1),1) histograms];
hcum=cumsum(h0,2);
medval=mean(hcum,2);
histp=hcum-medval;
histp=histp/max(hcum(:,end))*p.width;

if length(edges)>size(histograms,1)
    edges=edges(1:end-1)+(edges(2)-edges(1))/2;
end



% colors={'red','green','red','green'}
lutf=str2func(p.lut);
colors=lutf(numhist);
for k=1:numhist
    hvalh=[histp(:,k) histp(end:-1:1,k+1)];
    edgeh=[edges edges(end:-1:1)];
%     patch(hvalh(:), edgeh(:),colors{k})
    hp(k)=patch('XData',hvalh(:)+p.xposition, 'YData', edgeh(:),'FaceColor',colors(k,:));
%     hp(k).EdgeColor=hp.FaceColor;
    hp(k).FaceAlpha=p.FaceAlpha;
    hp(k).LineWidth=p.LineWidth;
    


end

    hvalh=[histp(:,1) histp(end:-1:1,end)];
    edgeh=[edges edges(end:-1:1)];
    hp(numhist+1)=patch('XData',hvalh(:)+p.xposition, 'YData', edgeh(:),'FaceColor','none');
    hp(numhist+1).LineWidth=p.LineWidth*2;



% statistics
datall=[];
for k=1:numhist
    if inputhist
        [mhist(k), nmed]=getstatisticshistogram(histograms(:,k),edges,p.statistics);
    else
        [mhist(k), nmed]=getstatisticsdat(dat{k},edges,p.statistics);
        datall=[datall; dat{k}];
    end
    yp(k)=mean(histp(nmed,k:k+1));
end
hold(ax,'on')
plot(ax,yp+p.xposition,mhist,'ko')

edgesall=repmat(edges,1,numhist);
if inputhist
    [mhistall, nmed,sall]=getstatisticshistogram(histograms(:),edgesall(:),p.statistics,p.quantile);
else
    [mhistall, nmed,sall]=getstatisticsdat(datall,edges,p.statistics,p.quantile);
end
errorbar(p.xposition,mhistall,sall(1),sall(2),'k','LineWidth',2,'CapSize',8);
errorbar(p.xposition,mhistall,0.00,'k','LineWidth',2,'CapSize',16)
end

function [mhist, nmed,s]=getstatisticshistogram(histograms,edges,statistics,quantile)
if nargin<4
    quantile= .95;
end
    switch statistics
        case 'mean'
            mhist=sum(histograms.*edges,1)./sum(histograms,1);
            s=sqrt(sum(histograms.*(edges-mhist).^2,1)./sum(histograms,1))*[-1 1];
            nmed=(find(mhist<=edges,1));
        case 'median'
            [esort,sind]=sort(edges);
               hc= cumsum(histograms(sind),1);
               nmed=find(hc(:)>=hc(end)/2,1);
               mhist=esort(nmed);
               ns=find(hc(:)>=hc(end)*(1-quantile),1);
               s(1)=esort(ns);
               ns=find(hc(:)>=hc(end)*(quantile),1);
               s(2)=esort(ns);
        otherwise
            disp('statistics should be mean or median')
    end

end

function [mhist, nmed,s]=getstatisticsdat(dat,edges,statistics,qv)
if nargin<4
    qv= .95;
end
    switch statistics
        case 'mean'
            mhist=mean(dat);
            s=std(dat)*[-1 1];
            
        case 'median'
            mhist=median(dat);
            s=quantile(dat,[1-qv qv]);
        otherwise
            disp('statistics should be mean or median')
    end
nmed=(find(mhist<=edges,1));
end


function [histograms, edges]=makehistogramsi(dat,binpositions)
if nargin>1 && ~isempty(binpositions)
    [histograms(:,1), edges]=histcounts(dat{1},binpositions);
else
    [histograms(:,1), edges]=histcounts(dat{1});
end
for k=2:length(dat)
    [histograms(:,k), edges]=histcounts(dat{k},edges);
end
edges=edges';
end


function po=parseinput(in)
p=inputParser;
addRequired(p,'binpositions');
addRequired(p,'histograms');
addOptional(p,'xposition',0);
addParameter(p,'Parent',[]);
addParameter(p,'lut','lines');
addParameter(p,'FaceAlpha',1);
addParameter(p,'LineWidth',1);
addParameter(p,'statistics','mean');
addParameter(p,'quantile',0.9);
addParameter(p,'width',1);
parse(p,in{:});
po=p.Results;
end