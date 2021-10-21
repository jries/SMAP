function modPoint = transModPoint(obj, modPoint, varargin)
% l stands for layer.
% Variation is now controled by SMAP as linkage error

p = inputParser;
p.addParameter('lPars',[]);
p.parse(varargin{:});
p = p.Results;
%% Transform the model
if isempty(p.lPars)
    lPars = obj.exportPars(1,'lPar');
else
    lPars = p.lPars;
end
if iscell(modPoint)
    for l = 1:obj.numOfLayer
        modPoint_.xnm = modPoint{l}.x; modPoint_.ynm = modPoint{l}.y;
        if isfield(modPoint{l}, 'z')
            modPoint_.znm = modPoint{l}.z;
        end
        modPoint_ = obj.locsHandler(modPoint_, lPars,0,'usedformalism', 'rotationMatrixRev','order_transform','RT');
        
        modPoint{l}.x = modPoint_.xnm; modPoint{l}.y = modPoint_.ynm;
        if isfield(modPoint{l}, 'z')
            modPoint{l}.z = modPoint_.znm;
        end
    end
else
    modPoint_.xnm = modPoint.x; modPoint_.ynm = modPoint.y;
    if isfield(modPoint, 'z')
        modPoint_.znm = modPoint.z;
    end
    modPoint_ = obj.locsHandler(modPoint_, lPars,0,'usedformalism', 'rotationMatrixRev','order_transform','RT');
    
    modPoint.x = modPoint_.xnm; modPoint.y = modPoint_.ynm;
    if isfield(modPoint, 'z')
        modPoint.z = modPoint_.znm;
    end
end
end