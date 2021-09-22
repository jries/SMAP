classdef summarize_NPCModelSelection<interfaces.DialogProcessor&interfaces.SEProcessor
    methods
        function obj=summarize_NPCModelSelection(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'se_viewer', 'se_siteroi', 'se_sitefov'};
            obj.showresults=true;
        end
        
        function out=run(obj,p)
            se = obj.locData.SE;
            sites = se.sites;

            list3 = getFieldAsVector(sites, 'annotation.list3.value');
            lUsed = list3==1;
%             siteOrder = 1:sum(lUsed);

           
            LLfit = [];
            AICc = [];
            normAICc = [];
            
            for m = 1:length(p.ID_LocMoFit)
                currentSym = p.symmetry(m);
                currentLocMoFitGUI = p.ID_LocMoFit(m);
                LLfit.(['sym',num2str(currentSym),'f']) = getFieldAsVector(sites(lUsed), ['evaluation.LocMoFitGUI_' num2str(currentLocMoFitGUI) '.fitInfo.LLfit']);
                AICc.(['sym',num2str(currentSym),'f']) = getFieldAsVector(sites(lUsed), ['evaluation.LocMoFitGUI_' num2str(currentLocMoFitGUI) '.fitInfo.AICc']);
                normAICc.(['sym',num2str(currentSym),'f']) = getFieldAsVector(sites(lUsed), ['evaluation.LocMoFitGUI_' num2str(currentLocMoFitGUI) '.fitInfo.normAICc']);
                LLExp.(['sym',num2str(currentSym),'f']) = getFieldAsVector(sites(lUsed), ['evaluation.LocMoFitGUI_' num2str(currentLocMoFitGUI) '.fitInfo.LLExp']);
                numOfLocsPerLayer.(['sym',num2str(currentSym),'f']) = getFieldAsVector(sites(lUsed), ['evaluation.LocMoFitGUI_' num2str(currentLocMoFitGUI) '.fitInfo.numOfLocsPerLayer']);
                LLExpMean.(['sym',num2str(currentSym),'f']) = getFieldAsVector(sites(lUsed), ['evaluation.LocMoFitGUI_' num2str(currentLocMoFitGUI) '.fitInfo.LLExpMean']);
                LLExpStd.(['sym',num2str(currentSym),'f']) = getFieldAsVector(sites(lUsed), ['evaluation.LocMoFitGUI_' num2str(currentLocMoFitGUI) '.fitInfo.LLExpStd']);
            end
            LLExp_gt.(['sym8f']) = getFieldAsVector(sites(lUsed), ['evaluation.gtLLExp.LLExp_gt']);
            LL_gt.(['sym8f']) = getFieldAsVector(sites(lUsed), ['evaluation.gtLLExp.LL_gt']);
            
            ax = obj.initaxis('Raw LL');
            ax_normAICc = obj.initaxis('normAICc');
            ax_normLL = obj.initaxis('normLL');
            axCmp = obj.initaxis('Comparison');
            axCmp_normAICc = obj.initaxis('Comparison_normAICc');
            axCmp_normLL = obj.initaxis('Comparison_normLL');
            ax_zscoreLL = obj.initaxis('zscoreLL');
            axCmp_zscoreLL = obj.initaxis('Comparison_zscoreLL');
            axCmp_LLExpTvsS = obj.initaxis('Comparison_LL_TvsS');
            axCmp_fit2gt = obj.initaxis('Comparison_fit2gt');
            axCmp_exp2data = obj.initaxis('Comparison_exp2data');
            axCmp_LL_gt2fit = obj.initaxis('Comparison_gt2fit');

%             colorID = [1 3 8 2 6];
            colorID = [2 8 1 6 3];
            %% LLfit: CDF
            axes(ax)
            hold on
%             palette= getPyPlot_cMap('tab10', 8,[],'"C:\Users\ries\AppData\Local\Programs\Python\Python37\python.exe"');
            for m = 1:length(p.ID_LocMoFit)
                currentSym = p.symmetry(m);
                currentLocMoFitGUI = p.ID_LocMoFit(m);
                LLfit_oneModel = LLfit.(['sym',num2str(currentSym),'f']);
                curve{m} = cdfplot(LLfit_oneModel);
                curve{m}.Color = myDiscreteLUT(colorID(m));
            end
            
            for m = length(p.ID_LocMoFit):-1:1
                currentSym = p.symmetry(m);
                legendText{m} = [num2str(currentSym) '-fold'];
            end
            xlabel(ax,'Maximum log-likelihood');
            ylabel(ax,'Cumulative probability');
            allLine = findobj(ax,'type','line');
            set(allLine,'linewidth',1.5)
            title(ax,[]);
            grid(ax, 'off')
            legend(legendText)
            hold off
            out = []; 
            
            %% normAICc: CDF
            axes(ax_normAICc)
            hold on
%             palette= getPyPlot_cMap('tab10', 8,[],'"C:\Users\ries\AppData\Local\Programs\Python\Python37\python.exe"');
            for m = 1:length(p.ID_LocMoFit)
                currentSym = p.symmetry(m);
                currentLocMoFitGUI = p.ID_LocMoFit(m);
                normAICc_oneModel = normAICc.(['sym',num2str(currentSym),'f']);
                curve_normAICc{m} = cdfplot(normAICc_oneModel);
                curve_normAICc{m}.Color = myDiscreteLUT(colorID(m));
            end
            
            xlabel(ax_normAICc,'Normalized AIC_c');
            ylabel(ax_normAICc,'Cumulative probability');
            allLine = findobj(ax_normAICc,'type','line');
            set(allLine,'linewidth',1.5)
            title(ax_normAICc,[]);
            grid(ax_normAICc, 'off')
            legend(legendText)
            hold off
            out = [];
            
            %% normLL: CDF
            axes(ax_normLL)
            hold on
%             palette= getPyPlot_cMap('tab10', 8,[],'"C:\Users\ries\AppData\Local\Programs\Python\Python37\python.exe"');
            for m = 1:length(p.ID_LocMoFit)
                currentSym = p.symmetry(m);
                currentLocMoFitGUI = p.ID_LocMoFit(m);
                normLL.(['sym',num2str(currentSym),'f']) = LLfit.(['sym',num2str(currentSym),'f'])-LLExp.(['sym',num2str(currentSym),'f']);
                normLL_oneModel = normLL.(['sym',num2str(currentSym),'f']);
                curve_normLL{m} = cdfplot(normLL_oneModel);
                curve_normLL{m}.Color = myDiscreteLUT(colorID(m));
            end
            
            xlabel(ax_normLL,'Normalized maximum log-likelihood');
            ylabel(ax_normLL,'Cumulative probability');
            allLine = findobj(ax_normLL,'type','line');
            set(allLine,'linewidth',1.5)
            title(ax_normLL,[]);
            grid(ax_normLL, 'off')
            legend(legendText)
            hold off
            out = [];
            %% zscore: CDF
            axes(ax_zscoreLL)
            hold on
%             palette= getPyPlot_cMap('tab10', 8,[],'"C:\Users\ries\AppData\Local\Programs\Python\Python37\python.exe"');
            for m = 1:length(p.ID_LocMoFit)
                currentSym = p.symmetry(m);
                currentLocMoFitGUI = p.ID_LocMoFit(m);
                zscoreLL.(['sym',num2str(currentSym),'f']) = (LLfit.(['sym',num2str(currentSym),'f'])-LLExpMean.(['sym',num2str(currentSym),'f']))./LLExpStd.(['sym',num2str(currentSym),'f']);
                zscoreLL_oneModel = zscoreLL.(['sym',num2str(currentSym),'f']);
                curve_zscoreLL{m} = cdfplot(zscoreLL_oneModel);
                curve_zscoreLL{m}.Color = myDiscreteLUT(colorID(m));
            end
            
            xlabel(ax_zscoreLL,'Z score');
            ylabel(ax_zscoreLL,'Cumulative probability');
            allLine = findobj(ax_zscoreLL,'type','line');
            set(allLine,'linewidth',1.5)
            title(ax_zscoreLL,[]);
            grid(ax_zscoreLL, 'off')
            legend(legendText)
            hold off
            out = [];
            
            %% LLFit: Comparison plot (6f vs 8f)
            pt = [LLfit.('sym8f');LLfit.('sym6f')]';
            Idx = rangesearch(pt,pt,0.1);
            count = cellfun(@length, Idx);
            scatter(axCmp, LLfit.('sym8f'), LLfit.('sym6f'),2, count, 'filled')
            hold(axCmp, 'on')
            plot(axCmp, [-17 -11],[-17 -11], '-k')
            hold(axCmp, 'off')
            xlabel(axCmp, 'Eight-fold symmetry')
            ylabel(axCmp, 'Six-fold symmetry')
            colorbar(axCmp)
            
            %% AICc: Comparison plot (6f vs 8f)
            pt = [normAICc.('sym8f');normAICc.('sym6f')]';
            Idx = rangesearch(pt,pt,0.1);
            count = cellfun(@length, Idx);
            scatter(axCmp_normAICc, normAICc.('sym8f'), normAICc.('sym6f'),2, count, 'filled')
            hold(axCmp_normAICc, 'on')
            plot(axCmp_normAICc, [25 31],[25 31], '-k')
            hold(axCmp_normAICc, 'off')
            xlabel(axCmp_normAICc, 'Eight-fold symmetry')
            ylabel(axCmp_normAICc, 'Six-fold symmetry')
            colorbar(axCmp_normAICc)
            
            %% NormLL: Comparison plot (6f vs 8f)
            pt = [normLL.('sym8f');normLL.('sym6f')]';
            Idx = rangesearch(pt,pt,0.1);
            count = cellfun(@length, Idx);
            scatter(axCmp_normLL, normLL.('sym8f'), normLL.('sym6f'),2, count, 'filled')
            hold(axCmp_normLL, 'on')
            plot(axCmp_normLL, [-0.2 0.2],[-0.2 0.2], '-k')
            hold(axCmp_normLL, 'off')
            xlabel(axCmp_normLL, 'Eight-fold symmetry')
            ylabel(axCmp_normLL, 'Six-fold symmetry')
            colorbar(axCmp_normLL)
            
            %% zscore: Comparison plot (6f vs 8f)
            pt = [normLL.('sym8f');normLL.('sym6f')]';
            Idx = rangesearch(pt,pt,0.1);
            count = cellfun(@length, Idx);
            scatter(axCmp_zscoreLL, zscoreLL.('sym8f'), zscoreLL.('sym6f'),2, count, 'filled')
            hold(axCmp_zscoreLL, 'on')
            plot(axCmp_zscoreLL, [-2 2],[-2 2], '-k')
            hold(axCmp_zscoreLL, 'off')
            xlabel(axCmp_zscoreLL, 'Eight-fold symmetry')
            ylabel(axCmp_zscoreLL, 'Six-fold symmetry')
            colorbar(axCmp_zscoreLL)
            
            %% LLExp: Comparison plot (theoretical vs simulated)
            scatter(axCmp_LLExpTvsS, LLExp.('sym8f'), LLExpMean.('sym8f'),5, 'filled')
            hold(axCmp_LLExpTvsS, 'on')
            plot(axCmp_LLExpTvsS, [-14.5 -11.5],[-14.5 -11.5], '-k')
            hold(axCmp_LLExpTvsS, 'off')
            xlabel(axCmp_LLExpTvsS, 'Theoretical LLFit')
            ylabel(axCmp_LLExpTvsS, 'Sampled LLFit')
            
            %% LLExp: Comparison plot (fit2gt)
            scatter(axCmp_fit2gt, LLExp.('sym8f'), LLExp_gt.('sym8f'),5, 'filled')
            hold(axCmp_fit2gt, 'on')
            plot(axCmp_fit2gt, [-14.5 -11.5],[-14.5 -11.5], '-k')
            hold(axCmp_fit2gt, 'off')
            xlabel(axCmp_fit2gt, 'Expected LL_{fit}')
            ylabel(axCmp_fit2gt, 'Expected LL_{gt}')
            
            %% Comparison plot (fit_exp2data)
            scatter(axCmp_exp2data, LLfit.('sym8f'), LLExp.('sym8f'),5, 'filled')
            hold(axCmp_exp2data, 'on')
            plot(axCmp_exp2data, [-14.5 -11.5],[-14.5 -11.5], '-k')
            hold(axCmp_exp2data, 'off')
            xlabel(axCmp_exp2data, 'LL_{fit}')
            ylabel(axCmp_exp2data, 'Expected LL_{fit}')
            
            %% LL: Comparison plot (gt2fit)
            scatter(axCmp_LL_gt2fit, LLfit.('sym8f'), LL_gt.('sym8f'),5, 'filled')
            hold(axCmp_LL_gt2fit, 'on')
            plot(axCmp_LL_gt2fit, [-14.5 -11.5],[-14.5 -11.5], '-k')
            hold(axCmp_LL_gt2fit, 'off')
            xlabel(axCmp_LL_gt2fit, 'LL_{fit}')
            ylabel(axCmp_LL_gt2fit, 'LL_{GT}')
            
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function recover_callback(a,b,obj)
    % get the pos from the backup
    obj.locData.loc.xnm = obj.locData.loc.xnmori;
    obj.locData.loc.ynm = obj.locData.loc.ynmori;
end

function rank_callback(a,b,obj)
    % get the pos from the backup
    fn = fieldnames(obj.locData.loc);
    a.String = fn';
    a.Enable = 'on';
end


function pard=guidef(obj)
pard.t1.object=struct('Style','text','String','ID (LocMoFitGUI)');
pard.t1.position=[1,1];
pard.t1.Width=2;
pard.ID_LocMoFit.object=struct('Style','edit','String','3 4 5 6 7');
pard.ID_LocMoFit.position=[1,3];
pard.ID_LocMoFit.Width=2;

pard.t2.object=struct('Style','text','String','Symmetry');
pard.t2.position=[2,1];
pard.t2.Width=2;
pard.symmetry.object=struct('Style','edit','String','6 7 8 9 10');
pard.symmetry.position=[2,3];
pard.symmetry.Width=2;
pard.plugininfo.type='ROI_Analyze';
pard.plugininfo.description='This plugin shows the summary of the model selection.';
end