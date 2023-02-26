classdef gtLLExp<interfaces.SEEvaluationProcessor
    %     calcuculate number of localizations with different dark times. Used
    %     for locsfromSE plugin (see counting)
    methods
        function obj=gtLLExp(varargin)
            obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function out=run(obj,p)
            %% Init
            out.LLExp_gt = [];
            %             layerson=obj.locData.parameters.sr_layerson;
            allParsArg = obj.site.evaluation.simulatesites.model.allParsArg;
            model = obj.site.evaluation.simulatesites.model.model;
            
            %% get and set up locs
            k = 1;
            
            fieldQ = {'locprecnm','locprecznm','xnm','znm','ynm','frame','xnmrot','ynmrot'};    % fields will be used
            fieldQNLayer = [fieldQ 'layer'];
            
            locs=obj.getLocs(fieldQ,'layer',k,'size','freeroi');
            locs.layer = ones(size(locs.(fieldQ{1})))*k;
            fieldExact = fieldnames(locs);
            lRm = ~ismember(fieldExact, fieldQNLayer);
            locs = rmfield(locs,fieldExact(lRm));
            fldElement = structfun(@(fld) length(fld),locs);
            fnLocs = fieldnames(locs);
            indFn2Rm = find(fldElement==0);
            locsOri = locs;
            locs = [];
            for l = length(fnLocs):-1:1
                if ~ismember(l,indFn2Rm)
                    locs.(fnLocs{l}) = locsOri.(fnLocs{l});
                end
            end
            locs.xnm = locs.xnmrot;
            locs.ynm = locs.ynmrot;
            
            %% set up the fitter
            fitter = LocMoFit('DataDim', 3);
            fitter.roiSize = obj.getPar('se_siteroi');
            
            fitter.locs = locs;
            
            fitter.addModel(model{1})
            fitter.getLocsInfo;
            fitter.advanceSetting.m91_background.value = 'density';

            %% Calculate
            fitter.setParArgBatch(allParsArg);
            fitter.setParArg('m1.lPar.variation','value',p.linkageError);

            % Expected LL at GT
            out.LLExp_gt = fitter.getELL([],[],5);

            % LL at GT
            out.LL_gt = fitter.loglikelihoodFun([],[],[])./sum(fitter.numOfLocsPerLayer);
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end

function pard=guidef
pard.text1.object=struct('Style','text','String','Linkage error');
pard.text1.position=[2,1];
pard.text1.Width=2;
pard.linkageError.object=struct('Style','edit','String','0');
pard.linkageError.position=[2,3];
pard.linkageError.Width=2;
pard.plugininfo.type='ROI_Evaluate';
pard.plugininfo.description='Evaluating the LLExp and LL at the ground truth.';
end