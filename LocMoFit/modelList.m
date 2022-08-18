function modelObj = modelList(modelName)
if nargin>0
switch modelName
case 'arc2D'
modelObj = arc2D;
case 'arc2D_arcLen'
modelObj = arc2D_arcLen;
case 'bucket2D'
modelObj = bucket2D;
case 'cspline3D_midPoint'
modelObj = cspline3D_midPoint;
case 'csplineClosedTube3D_midPoint'
modelObj = csplineClosedTube3D_midPoint;
case 'csplineTube3D_midPoint'
modelObj = csplineTube3D_midPoint;
case 'csplineTube3D_xyz'
modelObj = csplineTube3D_xyz;
case 'dualEllipse3D_avgR_discrete'
modelObj = dualEllipse3D_avgR_discrete;
case 'dualEllipse3D_discrete'
modelObj = dualEllipse3D_discrete;
case 'dualRing3D_discrete'
modelObj = dualRing3D_discrete;
case 'ellipse3D'
modelObj = ellipse3D;
case 'gaussianCluster2D'
modelObj = gaussianCluster2D;
case 'hemispheroid2D'
modelObj = hemispheroid2D;
case 'locsBG3D'
modelObj = locsBG3D;
case 'ring2D'
modelObj = ring2D;
case 'ring3D'
modelObj = ring3D;
case 'sphericalCap3D_surfaceArea'
modelObj = sphericalCap3D_surfaceArea;
case 'sphericalCap3Dp_surfaceArea'
modelObj = sphericalCap3Dp_surfaceArea;
case 'spheroid3Dp_surfaceArea'
modelObj = spheroid3Dp_surfaceArea;
case 'spheroidCap3Dp_surfaceArea'
modelObj = spheroidCap3Dp_surfaceArea;
case 'thickRing2D'
modelObj = thickRing2D;
end
else
modelObj = {'arc2D','arc2D_arcLen','bucket2D','cspline3D_midPoint','csplineClosedTube3D_midPoint','csplineTube3D_midPoint','csplineTube3D_xyz','dualEllipse3D_avgR_discrete','dualEllipse3D_discrete','dualRing3D_discrete','ellipse3D','gaussianCluster2D','hemispheroid2D','locsBG3D','ring2D','ring3D','sphericalCap3D_surfaceArea','sphericalCap3Dp_surfaceArea','spheroid3Dp_surfaceArea','spheroidCap3Dp_surfaceArea','thickRing2D'};
end
end
