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
            cutoffOneRing = 35;

            lUsed = getFieldAsVector(sites, 'annotation.use');
            siteOrder = 1:se.numberOfSites;

            usedSites = sites(lUsed);
            fn = fieldnames(obj.SE.processors.eval.children);
            if ismember('SMLMModelFitGUI_5', fn)
                whichSMLMModelFitGUI = '5';
            else
                whichSMLMModelFitGUI = '3';
            end
            fitInfo = getFieldAsVector(usedSites, ['evaluation.SMLMModelFitGUI_' whichSMLMModelFitGUI '.fitInfo']);
            lFailed = cellfun(@(x)strcmp(x.guiInfo,'Fit or plot failed.'), fitInfo);
            evalList = se.processors.eval.guihandles.modules.Data(:,2);
            indProcessor = find(strcmp('SMLMModelFitGUI',evalList));
            relativePosLastStep = 2;
            ID = getFieldAsVector(usedSites, 'ID');
            numOfLocsL1 = getFieldAsVectorInd(usedSites, ['evaluation.SMLMModelFitGUI_' whichSMLMModelFitGUI '.fitInfo.numOfLocsPerLayer'], 1);
            % [~,idxZ] = g.locData.SE.processors.eval.processors{indProcessor-2}.fitter.wherePar('pars.m1.lPar.z');
            % zS1 = getFieldAsVectorInd(usedSites, 'evaluation.SMLMModelFitGUI.allParsArg.value',idxZ);
            % 
            
            [~,idxRingDist] = se.processors.eval.processors{indProcessor}.fitter.wherePar('pars.m2.lPar.z');
            ringDistS1 = getFieldAsVectorInd(usedSites, 'evaluation.SMLMModelFitGUI.allParsArg.value',idxRingDist);

%             [cutoffOneRing, bin_edges] = getOneRingCutoff(ringDistS1);
            
            [~,idxZ] = se.processors.eval.processors{indProcessor+relativePosLastStep}.fitter.wherePar('pars.m1.lPar.z');
            z = getFieldAsVectorInd(usedSites, ['evaluation.SMLMModelFitGUI_' whichSMLMModelFitGUI '.allParsArg.value'],idxZ);

            [~,idxRingDist] = se.processors.eval.processors{indProcessor+relativePosLastStep}.fitter.wherePar('pars.m1.mPar.ringDistance');
            ringDist = getFieldAsVectorInd(usedSites, ['evaluation.SMLMModelFitGUI_' whichSMLMModelFitGUI '.allParsArg.value'],idxRingDist);

            [~,idxR] = se.processors.eval.processors{indProcessor+relativePosLastStep}.fitter.wherePar('pars.m1.mPar.radius');
            r = getFieldAsVectorInd(usedSites, ['evaluation.SMLMModelFitGUI_' whichSMLMModelFitGUI '.allParsArg.value'],idxR);

            [~,idxVar] = se.processors.eval.processors{indProcessor+relativePosLastStep}.fitter.wherePar('pars.m1.lPar.variation');
            var = getFieldAsVectorInd(usedSites, ['evaluation.SMLMModelFitGUI_' whichSMLMModelFitGUI '.allParsArg.value'],idxVar);

            [~,idxAzi] = se.processors.eval.processors{indProcessor+relativePosLastStep}.fitter.wherePar('pars.m1.mPar.azimuthalShift');
            azi = getFieldAsVectorInd(usedSites, ['evaluation.SMLMModelFitGUI_' whichSMLMModelFitGUI '.allParsArg.value'],idxAzi);
            azi = rem(rem(azi,45)+90, 45);
            
            [~,idxBG] = se.processors.eval.processors{indProcessor+relativePosLastStep}.fitter.wherePar('pars.m91.offset.weight');
            bg = getFieldAsVectorInd(usedSites, ['evaluation.SMLMModelFitGUI_' whichSMLMModelFitGUI '.allParsArg.value'],idxBG);
            
            numOfLocsLayer1 = getFieldAsVectorInd(usedSites, ['evaluation.SMLMModelFitGUI_' whichSMLMModelFitGUI '.fitInfo.numOfLocsPerLayer'],1);
            bgDensity = numOfLocsLayer1.*bg/(p.se_siteroi/1000)^3;
            
            lOneRing = ringDist<=cutoffOneRing;
            
            for k = find(lFailed)
                usedSites(k).annotation.list3.value = 5;
            end

            for k = find(lOneRing)'
                % 7: one-ring detected by the z correction
                % 8: one-ring detected by this plugin only
                % 9: one-ring detected by both
                if usedSites(k).annotation.list3.value == 7
                    usedSites(k).annotation.list3.value = 9;
                else
                    usedSites(k).annotation.list3.value = 8;
                end
            end
            list3 = getFieldAsVector(usedSites, 'annotation.list3.value');
            lGood = list3 == 1;
            
            %% Shift the athimuthal angle periodically 
            medAzi_0 = 0;
            medAzi = 90;
            
            azi_new = azi(lGood);
            while abs(medAzi_0-medAzi)>1e-6
                medAzi_0 = medAzi;
                medAzi = median(azi_new);
                azi_ub = medAzi+22.5;
                azi_lb = medAzi-22.5;
                azi_new(azi_new>azi_ub) = azi_new(azi_new>azi_ub)-45;
                azi_new(azi_new<azi_lb) = azi_new(azi_new<azi_lb)+45;
            end

            %%
            %% Ring separation
            axRS=obj.initaxis('Ring separation vs z');
            f = fit(z(lGood),ringDist(lGood),'poly1','Robust','LAR');
            h1 = plot(axRS, z(lGood),ringDist(lGood), ' ob');
            hold(axRS,'on')
            h2 = plot(axRS, z(~lGood),ringDist(~lGood), ' ok');
            h3 = plot(f, 'c');
            hold(axRS,'off')
            xlabel(axRS, 'z position (nm)')
            ylabel(axRS, 'Ring separation (nm)')
            legend([h1 h2 h3],{'Included','Excluded','Linear fit'})
            
            ax1 = obj.initaxis('Ring separation');
            par = ringDist(lGood);
            binWidth = 3;
            bin_Edge = floor(min(par)/binWidth)*binWidth:binWidth:ceil(max(par)/binWidth)*binWidth;
            histogram(ax1, par, bin_Edge)
            title(ax1, [sprintf('%.1f',mean(par)) '\pm' sprintf('%.1f', std(par))])
            xlabel(ax1, 'Distance (nm)')
            ylabel(ax1, 'Count')
            
            ax4 = obj.initaxis('Ring separation all');
            par = ringDist;
            binWidth = 3;
            bin_Edge = floor(min(par)/binWidth)*binWidth:binWidth:ceil(max(par)/binWidth)*binWidth;
            histogram(ax4, par, bin_Edge)
            title(ax4, [sprintf('%.1f',mean(par)) '\pm' sprintf('%.1f', std(par))])
            xlabel(ax4, 'Distance (nm)')
            ylabel(ax4, 'Count')
            
            ax2 = obj.initaxis('Ring twist');
            par = azi_new;
            binWidth = 4;
            bin_Edge = floor(min(par)/binWidth)*binWidth:binWidth:ceil(max(par)/binWidth)*binWidth;
            histogram(ax2, par, bin_Edge)
            title(ax2, [sprintf('%.1f',mean(par)) '\pm' sprintf('%.1f', std(par))])
            xlabel(ax2, ['Angle (' char(176) ')'])
            ylabel(ax2, 'Count')
            
            ax3 = obj.initaxis('Radius');
            par = r(lGood);
            binWidth = 2;
            bin_Edge = floor(min(par)/binWidth)*binWidth:binWidth:ceil(max(par)/binWidth)*binWidth;
            histogram(ax3, par, bin_Edge)
            title(ax3, [sprintf('%.1f',mean(par)) '\pm' sprintf('%.1f', std(par))])
            xlabel(ax3, 'Radius (nm)')
            ylabel(ax3, 'Count')
            
            ax4 = obj.initaxis('BGDensity');
            par = bgDensity(lGood);
            binWidth = 50;
            bin_Edge = floor(min(par)/binWidth)*binWidth:binWidth:ceil(max(par)/binWidth)*binWidth;
            histogram(ax4, par, bin_Edge)
            title(ax4, [sprintf('%.1f',mean(par)) '\pm' sprintf('%.1f', std(par))])
            xlabel(ax4, 'Density (um^-^3)')
            ylabel(ax4, 'Count')
            out = [];
            
            if p.showIntActPlot
                ax4 = obj.initaxis('Interactive'); 
                plotSElink(ax4,azi_new,numOfLocsL1,ID,se,' ob')
                xlabel(ax4, ['Twist angle (' char(176) ')'])
                ylabel(ax4, 'Number of Locs')
            end
            
            
            if p.showExampleMod
                fitter = se.processors.eval.processors{indProcessor+2}.fitter;
                fitter.setParArg('m1.mPar.ringDistance', 'value',50)
                fitter.setParArg('m1.mPar.azimuthalShift', 'value',8)
                fitter.setParArg('m1.mPar.radius', 'value',54)
                modPoints = fitter.model{1}.getPoint(fitter.exportPars(1,'mPar'));
                lLower = modPoints.z < 0;
                col = zeros([size(lLower,1) 3]);
                col(lLower,:) = repmat([0.9 0.9 0.9],sum(lLower),1);
                ax5 = obj.initaxis('Example model'); 
                scatter3(ax5, modPoints.x, modPoints.y, modPoints.z, 90, col,'filled','MarkerEdgeColor','k', 'LineWidth',1)
                axis(ax5,'equal')
                axis(ax5,'vis3d')
                set(ax5,'XLim',[-120 120])
                set(ax5,'YLim',[-120 120])
                view(0,45);
                grid(ax5, 'off')
                
            end
            
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

pard.showIntActPlot.object=struct('String','Display interactive','Style','checkbox');
pard.showIntActPlot.position=[2,1];
pard.showIntActPlot.Width=1.5;

pard.showExampleMod.object=struct('String','Display example model','Style','checkbox');
pard.showExampleMod.position=[3,1];
pard.showExampleMod.Width=1.5;

pard.plugininfo.type='ROI_Analyze';

end