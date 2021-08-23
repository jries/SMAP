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
            fitInfo = getFieldAsVector(usedSites, 'evaluation.LocMoFitGUI_3.fitInfo');
            evalList = se.processors.eval.guihandles.modules.Data(:,2);
            indProcessor = find(strcmp('LocMoFitGUI_3',evalList));
%             relativePosLastStep = 2;
            ID = getFieldAsVector(usedSites, 'ID');
            numOfLocsL1 = getFieldAsVectorInd(usedSites, 'evaluation.LocMoFitGUI_3.fitInfo.numOfLocsPerLayer', 1);
            % [~,idxZ] = g.locData.SE.processors.eval.processors{indProcessor-2}.fitter.wherePar('pars.m1.lPar.z');
            % zS1 = getFieldAsVectorInd(usedSites, 'evaluation.LocMoFitGUI.allParsArg.value',idxZ);
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
%                 oneGrpId = regexp(oriName, '\_(cond\d+\_\D+\d+)\_','tokens');
                oneGrpId = regexp(oriName, '_(\D+\d+)_','tokens');
                allGrpId(k) = oneGrpId;
            end
            allGrpId = repmat(allGrpId,[200 1]);
            allGrpId = [allGrpId{:}]';
            [~,~,stats] = errorPlot(obj, settings, parStack, allGrpId);
            
%             [~,idxVar] = se.processors.eval.processors{indProcessor+relativePosLastStep}.fitter.wherePar('pars.m1.lPar.variation');
%             var = getFieldAsVectorInd(usedSites, 'evaluation.LocMoFitGUI_3.allParsArg.value',idxVar);

            colID = [7 2];
            if p.separationHistogram
                axHist = obj.initaxis('Separation (LE)');
                allGrp2plot = {'cond2_LE0','cond2_LE4'};
                hold(axHist,'on')
                lParameter = strcmp({settings.tabTitle}, 'Separation');
                error = parStack(lParameter).fit-parStack(lParameter).gt;
                for l = 1:length(allGrp2plot)
                    lGrp = strcmp(allGrp2plot{l}, allGrpId);
                    histogram(axHist,error(lGrp)+50,-10:2.5:85, 'FaceColor',myDiscreteLUT(colID(l)))
                end
                hold(axHist,'off')
                xlabel(axHist,'Ring separation (nm)')
                ylabel(axHist,'Counts')
                axHist.LineWidth = 1.5;
                
                axScatter = obj.initaxis('Interactive');
                plotSElink(axScatter,error(lGrp)+50,1:sum(lGrp), ID(lGrp), se, ' ob')
            end
            
            try
                clipboard('copy', stats)
            catch
            end
            
            out = [];
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end
function settings = savedSettings()
    settings(1).parID = 'pars.m1.lPar.x';
    settings(1).tabTitle = 'Position \it x_0';
    settings(1).yLabUnit = 'nm';
    settings(1).bound = 10;

%     settings(end+1).parID = 'pars.m1.lPar.y';
%     settings(end).tabTitle = 'Y';
%     settings(end).yLabUnit = 'nm';
%     settings(end).bound = 20;

    settings(end+1).parID = 'pars.m1.lPar.z';
    settings(end).tabTitle = 'Position \it z_0';
    settings(end).yLabUnit = 'nm';
    settings(end).bound = 14;

    settings(end+1).parID = 'pars.m1.lPar.xrot';
    settings(end).tabTitle = 'Rotation {\alpha}';
    settings(end).yLabUnit = '\circ';
    settings(end).bound = 20;

    settings(end+1).parID = 'pars.m1.mPar.ringDistance';
    settings(end).tabTitle = 'Separation \it s';
    settings(end).yLabUnit = 'nm';
    settings(end).bound = 30;

    settings(end+1).parID = 'pars.m1.mPar.radius';
    settings(end).tabTitle = 'Radius \it r';
    settings(end).yLabUnit = 'nm';
    settings(end).bound = 5;

    settings(end+1).parID = 'pars.m1.mPar.azimuthalShift';
    settings(end).tabTitle = 'Twist {\theta}';
    settings(end).processFit = @(x) centerAroundGT(x, 8.8);
    settings(end).yLabUnit = '\circ';
    settings(end).bound = 12;

    settings(end+1).parID = 'pars.m91.offset.weight';
    settings(end).tabTitle = 'Background \it w_b_g';
    settings(end).processError = @(x) x*100;
    settings(end).yLabUnit = '%';
    settings(end).bound = 15;
    
    settings(end+1).parID = 'pars.m1.lPar.variation';
    settings(end).tabTitle = 'Linkage error {\epsilon}';
    settings(end).gt_alternative = 'evaluation.simulatesites.linkageerror';
    settings(end).yLabUnit = 'nm';
    settings(end).bound = 3;
end

function [grpMean, grpStd, table] = errorPlot(obj, settings, parStack, grp)
    generalTable = [];
    for k = 1:length(settings)
        oneSetting = settings(k);
        if size(parStack(k).gt,2)>1
            parStack(k).gt = parStack(k).gt';
        end
        error = parStack(k).fit-parStack(k).gt;
        if ~isempty(oneSetting.processError)
            error = oneSetting.processError(error);
        end
        [grpMean, grpStd, grpID] = grpstats(error,grp, {'mean','std','gname'});
%         grpID = regexp(grpID, 'cond(\d+)_(\D+)(\d+)','tokens');
        grpID = regexp(grpID, '(\D+)(\d+)','tokens');
%         cond = cellfun(@(x)str2num(x{1}{1}), grpID);
        varName = cellfun(@(x)x{1}{1}, grpID, 'UniformOutput', false);
        var = cellfun(@(x)str2num(x{1}{2}), grpID);
        
        % check which of BG, RB, and LE
        [uniVarName,~,ic] = unique(varName);
        varNameCounts = accumarray(ic,1);
        [~,maxID] = max(varNameCounts);
        [~,minID] = min(varNameCounts);
        testVar = uniVarName{maxID};
        
        if strcmp(testVar, 'RB')
            var(ic==minID) = 2;
        end
        
        [var,varOrder] = sort(var);
        grpMean = grpMean(varOrder);
        grpStd = grpStd(varOrder);
%         cond = cond(varOrder);
        
        ax = obj.initaxis(oneSetting.tabTitle);
        hold(ax,'on')
        h = {};
%         col = {'#ff8000','#008080','#0804a4' };
        col = '#000000';
%         conditionLabels = {'PALM','STORM','PAINT'};
        switch testVar
            case 'bg'
                var = var*10;
                xlab = 'Background (%)';
            case 'LE'
                var = ((1-var*0.1)-0.3)*100;
                xlab = 'Labeling efficiency (%)';
            case 'RB'
                xlab = 'Re-blinks';
            otherwise
        end
        
        if nargout>2
%             sampleLabels = strcat(testVar, '_', string(var) ,'-', conditionLabels(cond)');
            sampleLabels = strcat(testVar, '_', string(var));
            stats = strcat(strtrim(string(num2str(grpMean, '%0.1f'))), " ", char(177), " ", strtrim(string(num2str(grpStd, '%0.1f'))));
            table = [sampleLabels stats repelem(string(oneSetting.tabTitle), length(stats))' repelem(string(oneSetting.yLabUnit), length(stats))'];
            generalTable = [generalTable;table];
        end
        
            h = errorshade(ax, var, grpMean, grpStd,...
                'Color',col,...
                'LineWidth',1.5,...
                'Marker','o',...
                'MarkerFaceColor',col,...
                'Shade_FaceAlpha',0.25,...
                'Shade_FaceColor','none',...
                'Shade_LineStyle','--',...
                'Shade_EdgeColor',col,...
                'Shade_LineWidth',1);
%             h{con} = errorbar(ax, var(cond==con), grpMean(cond==con), grpStd(cond==con));

        lineObj = findobj(ax,'type','line');
        patchObj = findobj(ax,'type','Patch');
        ax.Children = [lineObj; patchObj];
%         legend([h{:}], {'PALM','STORM','PAINT'},'Location','southwest','NumColumns',3)
        hold(ax,'off');
%         range_var = range(var);
%         ax.XLim = [min(var)-0.1*range_var max(var)+0.1*range_var];
%         max_var = max(grpMean+0.6*grpStd);
%         min_var = min(grpMean-0.6*grpStd);
%         range_var = max_var-min_var;
%         ax.YLim = [round(min_var-range_var) round(max_var+range_var)];
        ax.YLim = [-oneSetting.bound oneSetting.bound];
        title(ax, oneSetting.tabTitle)
        xlabel(ax, xlab)
        ylabel(ax, ['Error (' oneSetting.yLabUnit ')'])
        
    end
    table = sprintf('%s\t%s\t%s\t%s\n', generalTable');
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
    val = getFieldAsVectorInd(sites, 'evaluation.LocMoFitGUI_3.allParsArg.value',idx);
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

pard.separationHistogram.object=struct('String','Histogram (ringSep)','Style','checkbox');
pard.separationHistogram.position=[1,1];
pard.separationHistogram.Width=1.5;

pard.plugininfo.type='ROI_Analyze';

end