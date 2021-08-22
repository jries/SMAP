classdef summarizeLocprecTitration<interfaces.DialogProcessor&interfaces.SEProcessor
    methods
        function obj=summarizeLocprecTitration(varargin)        
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
            
            
            switch strtrim(p.xAxis.String(p.xAxis.Value,:))
                case 'photon'
                    ph = regexp(getFieldAsVector(obj.SE.files,'name'), 'L\d*P(\d)*B0R2L2','tokens');
                    ph = [ph{:}];
                    ph = str2double(string(ph));
                    grp = grpstats(obj.locData.grouploc.locprecnm, ph(obj.locData.grouploc.filenumber),'median');

                    filenumber = getFieldAsVector(obj.SE.sites, 'info.filenumber');
                    nSites = grpstats(ph(filenumber)',ph(filenumber)','numel');
                    grp = repelem(grp,nSites);
                case 'file'
                    filenumber = getFieldAsVector(obj.SE.sites, 'info.filenumber');
                    grp = filenumber;
            end
            [~,~,stats] = errorPlot(obj, settings, parStack, grp);
            
%             [~,idxVar] = se.processors.eval.processors{indProcessor+relativePosLastStep}.fitter.wherePar('pars.m1.lPar.variation');
%             var = getFieldAsVectorInd(usedSites, 'evaluation.LocMoFitGUI_3.allParsArg.value',idxVar);

            colID = [7 2];
            
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
    settings(1).tabTitle = 'X';
    settings(1).yLabUnit = 'nm';
    settings(1).oneSideBound = 20;

%     settings(end+1).parID = 'pars.m1.lPar.y';
%     settings(end).tabTitle = 'Y';
%     settings(end).yLabUnit = 'nm';
%     settings(end).oneSideBound = 20;

    settings(end+1).parID = 'pars.m1.lPar.z';
    settings(end).tabTitle = 'Z';
    settings(end).yLabUnit = 'nm';
    settings(end).oneSideBound = 20;

    settings(end+1).parID = 'pars.m1.lPar.xrot';
    settings(end).tabTitle = 'Angle \alpha';
    settings(end).yLabUnit = '\circ';
    settings(end).oneSideBound = 25;

    settings(end+1).parID = 'pars.m1.mPar.ringDistance';
    settings(end).tabTitle = 'Separation';
    settings(end).yLabUnit = 'nm';
    settings(end).oneSideBound = 20;

    settings(end+1).parID = 'pars.m1.mPar.radius';
    settings(end).tabTitle = 'Radius';
    settings(end).yLabUnit = 'nm';
    settings(end).oneSideBound = 20;

    settings(end+1).parID = 'pars.m1.mPar.azimuthalShift';
    settings(end).tabTitle = 'Twist';
    settings(end).processFit = @(x) centerAroundGT(x, 8.8);
    settings(end).yLabUnit = '\circ';
    settings(end).oneSideBound = 22.5;

    settings(end+1).parID = 'pars.m91.offset.weight';
    settings(end).tabTitle = 'Background';
    settings(end).processError = @(x) x*100;
    settings(end).yLabUnit = '%';
    settings(end).oneSideBound = 25;
    
    settings(end+1).parID = 'pars.m1.lPar.variation';
    settings(end).tabTitle = 'Linkage error';
    settings(end).gt_alternative = 'evaluation.simulatesites.linkageerror';
    settings(end).yLabUnit = 'nm';
    settings(end).oneSideBound = 10;
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
                
        % check which of BG, RB, and LE
        ax = obj.initaxis(oneSetting.tabTitle);
        hold(ax,'on')
        col = '#ff8000';
        
        locprenm = str2double(string(grpID));
        
        if nargout>2
            stats = strcat(strtrim(string(num2str(grpMean, '%0.1f'))), " ", char(177), " ", strtrim(string(num2str(grpStd, '%0.1f'))));
            table = [locprenm stats repelem(string(oneSetting.tabTitle), length(stats))' repelem(string(oneSetting.yLabUnit), length(stats))'];
            generalTable = [generalTable;table];
        end
        
        errorshade(ax, locprenm, grpMean, grpStd,...
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

        hold(ax,'off');
%         range_var = range(var);
%         ax.XLim = [min(var)-0.1*range_var max(var)+0.1*range_var];
%         max_var = max(grpMean+0.6*grpStd);
%         min_var = min(grpMean-0.6*grpStd);
%         range_var = max_var-min_var;
%         ax.YLim = [round(min_var-range_var) round(max_var+range_var)];
        ax.YLim = [-oneSetting.oneSideBound oneSetting.oneSideBound];
        title(ax, oneSetting.tabTitle)
        xlabel(ax, 'Median localization precision (nm)')
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

pard.xAxis.object=struct('Style','popupmenu','String','file|photon');
pard.xAxis.position=[1,1];
pard.xAxis.Width=1.5;

pard.plugininfo.type='ROI_Analyze';

end