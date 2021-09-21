function layerPoint = getLayerPoint(obj, modSamplingF)
modPoint = obj.getModPoint(modSamplingF);
for ch = 1:length(obj.allModelLayer)
    layerPoint{ch}.x = [];
    layerPoint{ch}.y = [];
    layerPoint{ch}.z = [];
    layerPoint{ch}.n = [];
    indModelOneCh = find(obj.modelLayer == obj.allModelLayer(ch));
    for k = 1:length(indModelOneCh)
        indOneModel = indModelOneCh(k);
        onePoint = modPoint{indOneModel};
        layerPoint{ch}.x = [layerPoint{ch}.x; onePoint.x];
        layerPoint{ch}.y = [layerPoint{ch}.y; onePoint.y];
        if obj.model{1}.dimension == 3
            layerPoint{ch}.z = [layerPoint{ch}.z; onePoint.z];
        end
        layerPoint{ch}.n = [layerPoint{ch}.n; onePoint.n];
    end
end
end