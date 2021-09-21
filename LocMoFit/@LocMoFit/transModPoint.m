function modPoint = transModPoint(obj, modPoint)
% l stands for layer.
% Variation is now controled by SMAP as linkage error

%% Transform the model
lPar = obj.exportPars(1,'lPar');
for l = 1:obj.numOfLayer
    modPoint_.xnm = modPoint{l}.x; modPoint_.ynm = modPoint{l}.y;
    if isfield(modPoint{l}, 'z')
        modPoint_.znm = modPoint{l}.z;
    end
    modPoint_ = obj.locsHandler(modPoint_, lPar,0,'usedformalism', 'rotationMatrixRev','order_transform','RT');
    
    modPoint{l}.x = modPoint_.xnm; modPoint{l}.y = modPoint_.ynm;
    if isfield(modPoint{l}, 'z')
        modPoint{l}.z = modPoint_.znm;
    end
end
end