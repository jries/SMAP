function modelPoint = getModPoint(obj, modSamplingF)
for k = obj.numOfModel:-1:1
    if k==1
        oneMPars = obj.exportPars(k,'mPar');
        if isempty(modSamplingF)
            modelPoint{k} = obj.model{k}.getPoint(oneMPars);
        else
            modelPoint{k} = obj.model{k}.getPoint(oneMPars,'factor',modSamplingF);
        end
        
    else
        oneMPars = obj.exportPars(k,'mPar');
        if isempty(modSamplingF)
            ref = obj.model{k}.getPoint(oneMPars);
        else
            ref = obj.model{k}.getPoint(oneMPars,'factor',modSamplingF);
        end
        
        % translate lPar to mPar
        oneLPars = obj.exportPars(k,'lPar');          % get lPars
%         fn = fieldnames(oneLPars);
%         lFn2rm = ismember(fn,{'xscale','yscale'});
%         fn = fn(~lFn2rm);
%         for l = 1:length(fn)
%             oneLPars.(fn{l}) = -oneLPars.(fn{l});
%         end
        
        pseudoLocs.xnm = ref.x;
        pseudoLocs.ynm = ref.y;
        if obj.model{1}.dimension == 3
            pseudoLocs.znm = ref.z;
        end
        newPseudoLocs = obj.locsHandler(pseudoLocs, oneLPars,[],'order_transform', 'RT','usedformalism','rotationMatrixRev');
        ref.x = newPseudoLocs.xnm;
        ref.y = newPseudoLocs.ynm;
        if obj.model{1}.dimension == 3
            ref.z = newPseudoLocs.znm;
        end
        modelPoint{k} = ref;
    end
end
end