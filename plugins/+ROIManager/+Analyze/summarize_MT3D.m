classdef summarize_MT3D<interfaces.DialogProcessor&interfaces.SEProcessor
    methods
        function obj=summarize_MT3D(varargin)        
            obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'se_viewer', 'se_siteroi', 'se_sitefov'};
            obj.showresults=true;
        end
        
        function out=run(obj,p)
            se = obj.locData.SE;
            sites = se.sites;
            
            lUsed = getFieldAsVector(sites, 'annotation.use');

            usedSites = sites(lUsed);
            fitInfo = getFieldAsVector(usedSites, 'evaluation.SMLMModelFitGUI_2.fitInfo');
            lFailed = cellfun(@(x)strcmp(x.guiInfo,'Fit or plot failed.'), fitInfo);
            evalList = se.processors.eval.guihandles.modules.Data(:,2);
            indProcessor = find(strcmp('SMLMModelFitGUI_2',evalList));
            
            ID = getFieldAsVector(usedSites, 'ID');
            
            [~,idxR] = se.processors.eval.processors{indProcessor}.fitter.wherePar('pars.m1.mPar.r');
            r = getFieldAsVectorInd(usedSites, 'evaluation.SMLMModelFitGUI_2.allParsArg.value',idxR);

%             [cutoffOneRing, bin_edges] = getOneRingCutoff(ringDistS1);

%             for k = find(lOneRing)'
%                 % 7: one-ring detected by the z correction
%                 % 8: one-ring detected by this plugin only
%                 % 9: one-ring detected by both
%                 if usedSites(k).annotation.list3.value == 7
%                     usedSites(k).annotation.list3.value = 9;
%                 else
%                     usedSites(k).annotation.list3.value = 8;
%                 end
%             end
%             list3 = getFieldAsVector(usedSites, 'annotation.list3.value');
%             lGood = list3 == 1;
            
            %% Shift the athimuthal angle periodically 


            %%
            %% Parameters
            ax1 = obj.initaxis('Radius');
            par = r;
            binWidth = 1;
            bin_Edge = floor(min(par)/binWidth)*binWidth:binWidth:ceil(max(par)/binWidth)*binWidth;
            histogram(ax1, par, bin_Edge)
            title(ax1, [sprintf('%.1f',mean(par)) '\pm' sprintf('%.1f', std(par))])
            xlabel(ax1, 'Radius (nm)')
            ylabel(ax1, 'Count')
            out = [];
            
            clipboard('copy', sprintf('%.1f',median(r)))
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