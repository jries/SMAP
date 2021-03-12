function edges = getHistogramEdge(value, binSize)
% YW added. To make sure the edges covers the entire range of the value.
    minVal = min(value(:),[],1);
    maxVal = max(value(:),[],1);
    minVal = floor(minVal/binSize)*binSize;
    maxVal = ceil(maxVal/binSize)*binSize;
    edges = minVal:binSize:maxVal;
end