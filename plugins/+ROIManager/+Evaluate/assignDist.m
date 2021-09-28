classdef assignDist<interfaces.SEEvaluationProcessor
    %     calcuculate number of localizations with different dark times. Used
    %     for locsfromSE plugin (see counting)
    methods
        function obj=assignDist(varargin)
            obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function out=run(obj,p)
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
            se = obj.locData.SE;
            locMoFitter = se.processors.eval.processors{3}.fitter;
            allParsArg = obj.site.evaluation.LocMoFitGUI_2.allParsArg;
            locMoFitter.setParArgBatch(allParsArg)
            locMoFitter.locs = locs;
            
            lPar = locMoFitter.exportPars(1, 'lPar');
            locs = locMoFitter.locsHandler(locs, lPar);
            zOffset = obj.site.evaluation.LocMoFitGUI_2.fitInfo.modelPar_internal{1}.zOffset;
            locs.znm = locs.znm+zOffset;
            [azimuth,elevation,d] = cart2sph(locs.xnm,locs.ynm,locs.znm);
            
            r = abs(obj.site.evaluation.LocMoFitGUI_2.fitInfo.derivedPars{1}.radius);
            dd = d-r;
            lShow = dd<50;
            
            thetaFit = locMoFitter.getVariable('m1.closeAngle');
            lb = deg2rad(-thetaFit);
            lAbove = elevation>lb-deg2rad(10);
            
            
            
            [elevation_sorted, ind] = sort(elevation(lShow&lAbove));
            
            dd_shown = dd(lShow&lAbove);
            
            if obj.getPar('se_display')
                ax = obj.setoutput('Residue');
                scatter(ax,d(lShow).*cos(elevation(lShow)), locs.znm(lShow), [], dd(lShow));
                axis(ax,'equal')
                axP = obj.setoutput('Polar');
                plot(axP,elevation(lShow),dd(lShow), ' .', 'MarkerSize', 5);
                hold(axP, 'on')
                medDD = movmedian(dd_shown(ind),deg2rad(30), 'SamplePoints', double(elevation_sorted)+(1:length(elevation_sorted))'*eps);
                plot(axP, elevation_sorted, medDD, '-', 'LineWidth',2)
                hold(axP, 'off')
            end
            
%             
%             fullRange = range([r+eps lb-eps]);
%             locs_znm = locs.znm(lShow);
%             dd_show = dd(lShow);
            
            out.dd = dd(lShow);
%             out.distance2ref = d;
            out.ele = elevation(lShow);
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