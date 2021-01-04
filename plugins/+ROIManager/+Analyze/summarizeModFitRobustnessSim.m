classdef summarizeModFitRobustnessSim<interfaces.DialogProcessor&interfaces.SEProcessor
    methods
        function obj=summarizeModFitRobustnessSim(varargin)        
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
            evalList = se.processors.eval.guihandles.modules.Data(:,2);
            indProcessor = find(strcmp('SMLMModelFitGUI_3',evalList));
%             relativePosLastStep = 2;
            ID = getFieldAsVector(usedSites, 'ID');
            numOfLocsL1 = getFieldAsVectorInd(usedSites, 'evaluation.SMLMModelFitGUI_3.fitInfo.numOfLocsPerLayer', 1);
            % [~,idxZ] = g.locData.SE.processors.eval.processors{indProcessor-2}.fitter.wherePar('pars.m1.lPar.z');
            % zS1 = getFieldAsVectorInd(usedSites, 'evaluation.SMLMModelFitGUI.allParsArg.value',idxZ);
            %             
            fitter = se.processors.eval.processors{indProcessor}.fitter;
            
            % Calculate settings
            settings = savedSettings;
             
            parStack = extractPar(sites, fitter, settings);

%             settings(4).parID = 'pars.m1.lPar.variation';
%             settings(4).tabTitle = 'Linkage error';
%             
%             settings(5).parID = 'pars.m1.lPar.variation';
%             settings(5).tabTitle = 'Linkage error';
%             
%             settings(6).parID = 'pars.m1.lPar.variation';
%             settings(6).tabTitle = 'Linkage error';
            
            allGrpId = {};
            for k = 1:length(obj.locData.files.file)
                [~,oriName] = fileparts(obj.locData.files.file(k).name);
                oneGrpId = regexp(oriName, '\_(cond\d+\_\D+\d+)\_','tokens');
                allGrpId(k) = oneGrpId;
            end
            allGrpId = repmat(allGrpId,[200 1]);
            allGrpId = [allGrpId{:}]';
            errorPlot(obj, settings, parStack, allGrpId)
            
%             [~,idxVar] = se.processors.eval.processors{indProcessor+relativePosLastStep}.fitter.wherePar('pars.m1.lPar.variation');
%             var = getFieldAsVectorInd(usedSites, 'evaluation.SMLMModelFitGUI_3.allParsArg.value',idxVar);

            

            
            clipboard('copy', sprintf('%.1f',medAzi))
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end
function settings = savedSettings()
    settings(1).parID = 'pars.m1.lPar.x';
    settings(1).tabTitle = 'X';
    settings(1).yLabUnit = 'nm';

    settings(2).parID = 'pars.m1.lPar.y';
    settings(2).tabTitle = 'Y';
    settings(2).yLabUnit = 'nm';

    settings(3).parID = 'pars.m1.lPar.z';
    settings(3).tabTitle = 'Z';
    settings(3).yLabUnit = 'nm';

    settings(4).parID = 'pars.m1.lPar.xrot';
    settings(4).tabTitle = 'X tilt';
    settings(4).yLabUnit = '\circ';

    settings(5).parID = 'pars.m1.mPar.ringDistance';
    settings(5).tabTitle = 'Separation';
    settings(5).yLabUnit = 'nm';

    settings(6).parID = 'pars.m1.mPar.radius';
    settings(6).tabTitle = 'Radius';
    settings(6).yLabUnit = 'nm';

    settings(7).parID = 'pars.m1.mPar.azimuthalShift';
    settings(7).tabTitle = 'Twist';
    settings(7).processFit = @(x) centerAroundGT(x, 8.8);
    settings(7).yLabUnit = '\circ';

    settings(8).parID = 'pars.m91.offset.weight';
    settings(8).tabTitle = 'Background';
    settings(8).processError = @(x) x*100;
    settings(8).yLabUnit = '%';
    
    settings(9).parID = 'pars.m1.lPar.variation';
    settings(9).tabTitle = 'Linkage error';
    settings(9).gt_alternative = 'evaluation.simulatesites.linkageerror';
    settings(9).yLabUnit = 'nm';
end

function [grpMean, grpStd] = errorPlot(obj, settings, parStack, grp)
    for k = 1:length(settings)
        oneSetting = settings(k);
        error = parStack(k).fit-parStack(k).gt;
        if ~isempty(oneSetting.processError)
            error = oneSetting.processError(error);
        end
        [grpMean, grpStd, grpID] = grpstats(error,grp, {'mean','std','gname'});
        grpID = regexp(grpID, 'cond(\d+)_(\D+)(\d+)','tokens');
        cond = cellfun(@(x)str2num(x{1}{1}), grpID);
        varName = cellfun(@(x)x{1}{2}, grpID, 'UniformOutput', false);
        var = cellfun(@(x)str2num(x{1}{3}), grpID);
        ax = obj.initaxis(oneSetting.tabTitle);
        hold(ax,'on')
        h = {};
        switch varName{1}
            case 'bg'
                var = var*10;
                xlab = 'Background (%)';
            otherwise
        end
        for con = 1:max(cond)
            h{con} = errorbar(ax, var(cond==con), grpMean(cond==con), grpStd(cond==con));
        end
        legend([h{:}], {'PALM','STORM','PAINT'},'Location','southwest','NumColumns',3)
        hold(ax,'off');
        range_var = range(var);
        ax.XLim = [min(var)-0.1*range_var max(var)+0.1*range_var];
        title(ax, oneSetting.tabTitle)
        xlabel(ax, xlab)
        ylabel(ax, ['Error (' oneSetting.yLabUnit ')'])
    end

end

function parStack = extractPar(sites, fitter, settings)
    parStack = [];
    for k = 1:length(settings)
        [val,val_gt] = extractOnePar(sites, fitter, settings(k).parID, settings(k).gt_alternative);
        parStack(k).fit = val;
        parStack(k).gt = val_gt;
        if ~isempty(settings(k).processFit)
            parStack(k).fit = settings(k).processFit(parStack(k).fit);
        end
    end
end

function [val,val_gt] = extractOnePar(sites, fitter, parID, gt_alternative)
    [~,idx] = fitter.wherePar(parID);
    val = getFieldAsVectorInd(sites, 'evaluation.SMLMModelFitGUI_3.allParsArg.value',idx);
    if isempty(gt_alternative)
        val_gt = getFieldAsVectorInd(sites, 'evaluation.simulatesites.model.allParsArg.value',idx);
    else
        val_gt = getFieldAsVector(sites, 'evaluation.simulatesites.linkageerror');
    end
end

function azi_new = centerAroundGT(azi, gt)
    azi = rem(rem(azi,45)+90, 45);
    medAzi_0 = gt;
    azi_new = azi;
    azi_ub = medAzi_0+22.5;
    azi_lb = medAzi_0-22.5;
    azi_new(azi_new>azi_ub) = azi_new(azi_new>azi_ub)-45;
    azi_new(azi_new<azi_lb) = azi_new(azi_new<azi_lb)+45;
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