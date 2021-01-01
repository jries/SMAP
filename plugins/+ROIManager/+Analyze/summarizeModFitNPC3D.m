classdef summarizeModFitNPC3D<interfaces.DialogProcessor&interfaces.SEProcessor
    methods
        function obj=summarizeModFitNPC3D(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'se_viewer', 'se_siteroi', 'se_sitefov'};
            obj.showresults=true;
        end
        
        function out=run(obj,p)
            se = obj.locData.SE;
            sites = se.sites;

            lUsed = getFieldAsVector(sites, 'annotation.use');
            siteOrder = 1:se.numberOfSites;

            usedSites = sites(lUsed);
            fitInfo = getFieldAsVector(usedSites, 'evaluation.SMLMModelFitGUI_3.fitInfo');
            lFailed = cellfun(@(x)strcmp(x.guiInfo,'Fit or plot failed.'), fitInfo);
            evalList = se.processors.eval.guihandles.modules.Data(:,2);
            indProcessor = find(strcmp('SMLMModelFitGUI',evalList));
            relativePosLastStep = 2;
            ID = getFieldAsVector(usedSites, 'ID');

            % [~,idxZ] = g.locData.SE.processors.eval.processors{indProcessor-2}.fitter.wherePar('pars.m1.lPar.z');
            % zS1 = getFieldAsVectorInd(usedSites, 'evaluation.SMLMModelFitGUI.allParsArg.value',idxZ);
            % 
            [~,idxRingDist] = se.processors.eval.processors{indProcessor}.fitter.wherePar('pars.m2.lPar.z');
            ringDistS1 = getFieldAsVectorInd(usedSites, 'evaluation.SMLMModelFitGUI.allParsArg.value',idxRingDist);

            [cutoffOneRing, bin_edges] = getOneRingCutoff(ringDistS1);
            
            [~,idxZ] = se.processors.eval.processors{indProcessor+relativePosLastStep}.fitter.wherePar('pars.m1.lPar.z');
            z = getFieldAsVectorInd(usedSites, 'evaluation.SMLMModelFitGUI_3.allParsArg.value',idxZ);

            [~,idxRingDist] = se.processors.eval.processors{indProcessor+relativePosLastStep}.fitter.wherePar('pars.m1.mPar.ringDistance');
            ringDist = getFieldAsVectorInd(usedSites, 'evaluation.SMLMModelFitGUI_3.allParsArg.value',idxRingDist);

            [~,idxR] = se.processors.eval.processors{indProcessor+relativePosLastStep}.fitter.wherePar('pars.m1.mPar.radius');
            r = getFieldAsVectorInd(usedSites, 'evaluation.SMLMModelFitGUI_3.allParsArg.value',idxR);

            [~,idxVar] = se.processors.eval.processors{indProcessor+relativePosLastStep}.fitter.wherePar('pars.m1.lPar.variation');
            var = getFieldAsVectorInd(usedSites, 'evaluation.SMLMModelFitGUI_3.allParsArg.value',idxVar);

            [~,idxAzi] = se.processors.eval.processors{indProcessor+relativePosLastStep}.fitter.wherePar('pars.m1.mPar.azimuthalShift');
            azi = getFieldAsVectorInd(usedSites, 'evaluation.SMLMModelFitGUI_3.allParsArg.value',idxAzi);
            azi = rem(rem(azi,45)+90, 45);
            
            %% Shift the athimuthal angle periodically 
            medAzi_0 = 0;
            medAzi = 90;
            lOneRing = ringDist<=cutoffOneRing;
            azi_new = azi(~(lOneRing|lFailed'));
            while abs(medAzi_0-medAzi)>1e-6
                medAzi_0 = medAzi;
                medAzi = median(azi_new);
                azi_ub = medAzi+22.5;
                azi_lb = medAzi-22.5;
                azi_new(azi_new>azi_ub) = azi_new(azi_new>azi_ub)-45;
                azi_new(azi_new<azi_lb) = azi_new(azi_new<azi_lb)+45;
            end

            %%
            for k = 1:length(usedSites)
                usedSites(k).annotation.list3.value = 1;
            end

            for k = find(lFailed)
                usedSites(k).annotation.list3.value = 9;
            end

            for k = find(lOneRing)'
                usedSites(k).annotation.list3.value = 8;
            end

            ax1 = obj.initaxis('Ring separation');
            par = ringDist(~lOneRing);
            binWidth = 3;
            bin_Edge = floor(min(par)/binWidth)*binWidth:binWidth:ceil(max(par)/binWidth)*binWidth;
            histogram(ax1, par, bin_Edge)
            title(ax1, [sprintf('%.1f',mean(par)) '\pm' sprintf('%.1f', std(par))])
            xlabel(ax1, 'Distance (nm)')
            ylabel(ax1, 'Count')
            
            ax2 = obj.initaxis('Ring twist');
            par = azi_new;
            binWidth = 4;
            bin_Edge = floor(min(par)/binWidth)*binWidth:binWidth:ceil(max(par)/binWidth)*binWidth;
            histogram(ax2, par, bin_Edge)
            title(ax2, [sprintf('%.1f',mean(par)) '\pm' sprintf('%.1f', std(par))])
            xlabel(ax2, ['Angle (' char(176) ')'])
            ylabel(ax2, 'Count')
            
            ax3 = obj.initaxis('Radius');
            par = r(~lOneRing);
            binWidth = 2;
            bin_Edge = floor(min(par)/binWidth)*binWidth:binWidth:ceil(max(par)/binWidth)*binWidth;
            histogram(ax3, par, bin_Edge)
            title(ax3, [sprintf('%.1f',mean(par)) '\pm' sprintf('%.1f', std(par))])
            xlabel(ax3, 'Radius (nm)')
            ylabel(ax3, 'Count')
            out = [];
            
            clipboard('copy', sprintf('%.1f',medAzi))
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