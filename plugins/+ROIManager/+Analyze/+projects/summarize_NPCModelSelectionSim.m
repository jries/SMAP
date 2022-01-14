classdef summarize_NPCModelSelectionSim<interfaces.DialogProcessor&interfaces.SEProcessor
    methods
        function obj=summarize_NPCModelSelectionSim(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'se_viewer', 'se_siteroi', 'se_sitefov'};
            obj.showresults=true;
        end
        
        function out=run(obj,p)
            se = obj.locData.SE;
            sites = se.sites;

            lUsed = getFieldAsVector(sites, 'annotation.use');
            siteOrder = 1:se.numberOfSites;

            [~,fName] = fileparts({obj.locData.files.file.name});
            sym = regexp(fName,'_(\d+)f_', 'tokens');
            sym = [sym{:}];
            sym = [sym{:}];
            sym = str2double(sym);
            sym = repelem(sym, 200);
%             for k = se.numberOfSites:-1:1
% %                 sym(k) = se.sites(k).evaluation.simulatesites.model.model{1}.modelObj.internalSettings.cornerNum;
%                 sym(k) = se.sites(k).evaluation.simulatesites.model.model{1}.modelObj.internalSettings.cornerNum;
%             end

            normAICc = [];
            
            for m = 1:5
                normAICc.(['sym',num2str(m+5),'f']) = getFieldAsVector(sites, ['evaluation.LocMoFitGUI_' num2str(m+2) '.fitInfo.normAICc']);
            end
            
            
            allSym = unique(sym);
            colorID = [2 8 1 6 3];
            for d = 1:length(allSym)
                ax = obj.initaxis(['data_sym',num2str(d+5),'f']);
                lOneSym = sym == allSym(d);
                axes(ax)
                hold on
                for m = 1:5
                    normAICc_oneModel = normAICc.(['sym',num2str(m+5),'f'])(lOneSym);
                    cdfplot(normAICc_oneModel)
                    curve{m} = cdfplot(normAICc_oneModel);
                    curve{m}.Color = myDiscreteLUT(colorID(m));
                end
                xlabel(ax,'Normalized AICC');
                ylabel(ax,'Cumulative probability');
                allLine = findobj(ax,'type','line');
                set(allLine,'linewidth',1.5)
                title(ax,[]);
                grid(ax, 'off')
                if d == 1
                    legend([curve{:}], {'6-fold','7-fold','8-fold','9-fold','10-fold'})
                end
%                 hold off
%                 out = [];    
                hold off
            end
            out = [];
            %% Comparison plot (6f vs 8f)
            axCmp = obj.initaxis('Comparison');
            pt = [normAICc.('sym8f')(sym == 8);normAICc.('sym6f')(sym == 8)]';
            Idx = rangesearch(pt,pt,0.1);
            count = cellfun(@length, Idx);
            scatter(axCmp, normAICc.('sym8f')(sym == 8), normAICc.('sym6f')(sym == 8),2, count, 'filled')
            hold(axCmp, 'on')
            plot(axCmp, [20 33],[20 33], '-k')
            hold(axCmp, 'off')
            xlabel(axCmp, 'Eight-fold symmetry')
            ylabel(axCmp, 'Six-fold symmetry')
            colorbar(axCmp)
            %% Mix (6f and 8f)
            axMix = obj.initaxis('Mix');
            sz = 4;
            hs8 = scatter(axMix, normAICc.('sym8f')(sym==8), normAICc.('sym6f')(sym==8),sz, myDiscreteLUT(colorID(3), 'type', 'rgb'), 'filled');
            hold(axMix, 'on')
            hs6 = scatter(axMix, normAICc.('sym8f')(sym==6), normAICc.('sym6f')(sym==6),sz, myDiscreteLUT(colorID(1), 'type', 'rgb'), 'filled');
            set([hs8 hs6],'MarkerFaceAlpha',0.3)
            plot(axMix, [20 33],[20 33], '-k')
            hold(axMix, 'off')
            xlabel(axMix, 'Eight-fold symmetry')
            ylabel(axMix, 'Six-fold symmetry')
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

pard.showIntActPlot.object=struct('String','Display interactive','Style','checkbox');
pard.showIntActPlot.position=[2,1];
pard.showIntActPlot.Width=1.5;

pard.showExampleMod.object=struct('String','Display example model','Style','checkbox');
pard.showExampleMod.position=[3,1];
pard.showExampleMod.Width=1.5;

pard.plugininfo.type='ROI_Analyze';

end