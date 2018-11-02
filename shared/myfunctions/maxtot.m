function [maxv,minv]=maxtot(in)

if iscell(in)
    maxv=-inf;
    minv=inf;
    for k=1:length(in)
        maxv=max(maxv,max(in{k}(:)));
        minv=min(minv,min(in{k}(:)));
    end
else
    maxv=(max(in(:)));
    minv=(min(in(:)));
end
        
        