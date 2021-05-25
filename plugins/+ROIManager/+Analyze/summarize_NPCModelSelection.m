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

            lUsed = getFieldAsVector(sites, 'annotation.use');
            siteOrder = 1:se.numberOfSites;

           
            LLfit = [];
            
            for m = 1:5
                LLfit.(['sym',num2str(m+5),'f']) = getFieldAsVector(sites, ['evaluation.SMLMModelFitGUI_' num2str(m+2) '.fitInfo.LLfit']);
            end
            
            ax = obj.initaxis('Raw LL');
            axes(ax)
            colorID = [1 3 8 2 6];
            hold on
            
%             palette= getPyPlot_cMap('tab10', 8,[],'"C:\Users\ries\AppData\Local\Programs\Python\Python37\python.exe"');
            for m = 1:5
                LLfit_oneModel = LLfit.(['sym',num2str(m+5),'f']);
                curve{m} = cdfplot(LLfit_oneModel);
                curve{m}.Color = myDiscreteLUT(colorID(m));
            end
            xlabel(ax,'Log-likelihood');
            ylabel(ax,'Cumulative probability');
            allLine = findobj(ax,'type','line');
            set(allLine,'linewidth',1.5)
            title(ax,[]);
            grid(ax, 'off')
            legend({'6-fold','7-fold','8-fold','9-fold','10-fold'})
            hold off
            out = [];    
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