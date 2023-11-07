function modImage = getModImage(obj, modelID, varargin)
% parse varargin
p = inputParser;
p.addParameter('mPars', obj.exportPars(modelID, 'mPars'));
p.addParameter('pixelSize', obj.model{modelID}.pixelSize);
p.addParameter('roiSize', obj.roiSize);

parse(p,varargin{:});
results = p.Results;
pixelSize = results.pixelSize;
roiSize = results.roiSize;

modImage = obj.model{modelID}.getImage(mPars, 'pixelSize', pixelSize, 'roiSize', roiSize);

if 1/pixelSize~=1
    if obj.dimension == 2
        modImage = imresize(modImage, 1/(pixelSize/obj.pixelSize));
    else
        modImage = imresize3(modImage, 1/(pixelSize/obj.pixelSize));
    end
end
end