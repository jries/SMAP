function stop = optimOutputPar(obj,x,optimValues,state)
    stop = false;
    switch state
        case 'init'
            obj.fitInfo.optimHistory = obj.parsInit.init(~obj.allParsArg.fix)';
            obj.fitInfo.optimHistory(end+1,:) = x;
        otherwise
            obj.fitInfo.optimHistory(end+1,:) = x;
    end
end