function stop = optimoOutputPar(obj,x,optimValues,state,varargin)
    stop = false;
    
    lFix = obj.allParsArg.fix;
    values_ori = obj.allParsArg.value(~lFix);
    x = xtransform(x,varargin{1});
    obj.allParsArg.value(~lFix) = x;
    
    [~,modCoord] = obj.plot(obj.locs,'plotType','point','modelSamplingFactor',0.3,'doNotPlot',true);
    lPars = obj.exportPars(1,'lPar');
    lPars.variation = 0;
    locsViz = obj.locsHandler(obj.locs,lPars,1);
    
    cla(gca)
    obj.rotCoordNMkImg(gca, modCoord, locsViz, [0 0], 10, 'Model', inf, {'red hot'});
%     obj.optFrames(size())
    oneFrame = getframe(gca);
    oneFrame = oneFrame.cdata;
    oneFrame = imresize(oneFrame,size(oneFrame,1:2)/4);
    if strcmp(state,'init')
        obj.optFrames = uint8(zeros([size(oneFrame) 5000]));
    end
    obj.optFrames(:,:,:,optimValues.iteration+1) = oneFrame;
    obj.allParsArg.value(~lFix) = values_ori;
    
end


function xtrans = xtransform(x,params)
% converts unconstrained variables into their original domains

xtrans = zeros(params.xsize);
% k allows some variables to be fixed, thus dropped from the
% optimization.
k=1;
for i = 1:params.n
  switch params.BoundClass(i)
    case 1
      % lower bound only
      xtrans(i) = params.LB(i) + x(k).^2;
      
      k=k+1;
    case 2
      % upper bound only
      xtrans(i) = params.UB(i) - x(k).^2;
      
      k=k+1;
    case 3
      % lower and upper bounds
      xtrans(i) = (sin(x(k))+1)/2;
      xtrans(i) = xtrans(i)*(params.UB(i) - params.LB(i)) + params.LB(i);
      % just in case of any floating point problems
      xtrans(i) = max(params.LB(i),min(params.UB(i),xtrans(i)));
      
      k=k+1;
    case 4
      % fixed variable, bounds are equal, set it at either bound
      xtrans(i) = params.LB(i);
    case 0
      % unconstrained variable.
      xtrans(i) = x(k);
      
      k=k+1;
  end
end

end % sub function xtransform end
