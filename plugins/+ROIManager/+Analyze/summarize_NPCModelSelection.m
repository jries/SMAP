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
            
            for m = 1:5
                LLfit.(['sym',num2str(m+5),'f']) = getFieldAsVector(sites(lUsed), ['evaluation.LocMoFitGUI_' num2str(m+2) '.fitInfo.LLfit']);
                AICc.(['sym',num2str(m+5),'f']) = getFieldAsVector(sites(lUsed), ['evaluation.LocMoFitGUI_' num2str(m+2) '.fitInfo.AICc']);
                normAICc.(['sym',num2str(m+5),'f']) = getFieldAsVector(sites(lUsed), ['evaluation.LocMoFitGUI_' num2str(m+2) '.fitInfo.normAICc']);
            end
            
            ax = obj.initaxis('Raw LL');
            ax_normAICc = obj.initaxis('normAICc');
            axCmp = obj.initaxis('Comparison');
            axCmp_normAICc = obj.initaxis('Comparison_normAICc');
%             colorID = [1 3 8 2 6];
            colorID = [2 8 1 6 3];
            %% LLfit: CDF
            axes(ax)
            hold on
%             palette= getPyPlot_cMap('tab10', 8,[],'"C:\Users\ries\AppData\Local\Programs\Python\Python37\python.exe"');
            for m = 1:5
                LLfit_oneModel = LLfit.(['sym',num2str(m+5),'f']);
                curve{m} = cdfplot(LLfit_oneModel);
                curve{m}.Color = myDiscreteLUT(colorID(m));
            end
            
           
            xlabel(ax,'Maximum Log-likelihood');
            ylabel(ax,'Cumulative probability');
            allLine = findobj(ax,'type','line');
            set(allLine,'linewidth',1.5)
            title(ax,[]);
            grid(ax, 'off')
            legend({'6-fold','7-fold','8-fold','9-fold','10-fold'})
            hold off
            out = []; 
            
            %% normAICc: CDF
            axes(ax_normAICc)
            hold on
%             palette= getPyPlot_cMap('tab10', 8,[],'"C:\Users\ries\AppData\Local\Programs\Python\Python37\python.exe"');
            for m = 1:5
                normAICc_oneModel = normAICc.(['sym',num2str(m+5),'f']);
                curve_normAICc{m} = cdfplot(normAICc_oneModel);
                curve_normAICc{m}.Color = myDiscreteLUT(colorID(m));
            end
            
            xlabel(ax_normAICc,'Normalized AIC_c');
            ylabel(ax_normAICc,'Cumulative probability');
            allLine = findobj(ax_normAICc,'type','line');
            set(allLine,'linewidth',1.5)
            title(ax_normAICc,[]);
            grid(ax_normAICc, 'off')
            legend({'6-fold','7-fold','8-fold','9-fold','10-fold'})
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

pard.plugininfo.type='ROI_Analyze';

end