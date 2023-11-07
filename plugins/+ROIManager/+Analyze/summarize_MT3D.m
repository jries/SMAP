classdef summarize_MT3D<interfaces.DialogProcessor&interfaces.SEProcessor
    methods
        function obj=summarize_MT3D(varargin)
            obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'se_viewer', 'se_siteroi', 'se_sitefov'};
            obj.showresults=true;
            
        end
        
        function makeGui(obj,varargin)
            makeGui@interfaces.DialogProcessor(obj); %make the main GUI from the guidef definitions
            %Settings
            obj.createGrpTable;
            comp_callback(obj.guihandles.comp,[],obj);
        end
        
        function out=run(obj,p)
            se = obj.locData.SE;
            sites = se.sites;
            
            lUsed = getFieldAsVector(sites, 'annotation.use');
            usedSites = sites(lUsed);
            fitInfo = getFieldAsVector(usedSites, 'evaluation.LocMoFitGUI_2.fitInfo');
            lFailed = cellfun(@(x)strcmp(x.guiInfo,'Fit or plot failed.'), fitInfo);
            lGood = getFieldAsVector(sites,'annotation.list3.value')==1;
            evalList = se.processors.eval.guihandles.modules.Data(:,2);
            indProcessor = find(strcmp('LocMoFitGUI_2',evalList));
            
            ID = getFieldAsVector(usedSites, 'ID');
            
            [~,idxR] = se.processors.eval.processors{indProcessor}.fitter.wherePar('pars.m1.mPar.r');
            r = getFieldAsVectorInd(usedSites, 'evaluation.LocMoFitGUI_2.allParsArg.value',idxR);
            
            [~,idxVar] = se.processors.eval.processors{indProcessor}.fitter.wherePar('pars.m1.lPar.variation');
            var = getFieldAsVectorInd(usedSites, 'evaluation.LocMoFitGUI_2.allParsArg.value',idxVar);
            
            [~,idxZMid] = se.processors.eval.processors{indProcessor}.fitter.wherePar('pars.m1.mPar.zMid');
            zMid = getFieldAsVectorInd(usedSites, 'evaluation.LocMoFitGUI_2.allParsArg.value',idxZMid);
                        
            if p.comp
                fileNumber = getFieldAsVector(usedSites, 'info.filenumber');
                grpLabel = p.grpTable.Data(:,2);
            end
            
            if ~isempty(usedSites(1).evaluation.LocMoFitGUI_2.fitInfo.derivedPars)&&isfield(usedSites(1).evaluation.LocMoFitGUI_2.fitInfo.derivedPars{1},'avgCurvature')
                for k = sum(lUsed):-1:1
                    curvature(k) = usedSites(k).evaluation.LocMoFitGUI_2.fitInfo.derivedPars{1}.avgCurvature;
                end
            end
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
            if p.comp
                hold(ax1, 'on')
                for k = 1:length(grpLabel)
                    if ~isempty(grpLabel{k})
                        histogram(ax1, par(fileNumber == k), bin_Edge);
                        finalLabel{k} = [grpLabel{k} ';' sprintf('%.1f',mean(par(fileNumber == k))) '\pm' sprintf('%.1f', std(par(fileNumber == k)))];
                    end
                end
                hold(ax1, 'off')
                legend(ax1,finalLabel(~cellfun(@isempty, finalLabel)))
            else
                histogram(ax1, par(lGood), bin_Edge);
                title(ax1, [sprintf('%.1f',mean(par(lGood))) '\pm' sprintf('%.1f', std(par(lGood)))])
            end
                xlabel(ax1, 'Radius (nm)')
                ylabel(ax1, 'Count')
                out = [];
                
                clipboard('copy', sprintf('%.1f',median(r)))
            
            ax2 = obj.initaxis('Variation vs r');
            par1 = r;
            par2 = var;
            
            if p.comp
                hold(ax2, 'on')
                for k = 1:length(grpLabel)
                    plot(ax2, par1(fileNumber == k), par2(fileNumber == k), ' o');
                end
                hold(ax2, 'off')
            else
                plot(ax2, par1, par2, ' o');
            end
            
            xlabel(ax2, 'Radius (nm)')
            ylabel(ax2, 'Variation (nm)')

            ax3 = obj.initaxis('z vs r');
            par1 = zMid;
            par2 = r;
            if p.comp
                hold(ax3, 'on')
                for k = 1:length(grpLabel)
                    plot(ax3, par1(fileNumber == k), par2(fileNumber == k), ' o');
                end
                hold(ax3, 'off')
            else
                plot(ax3, par1, par2, ' o');
            end

            xlabel(ax3, 'z position (nm)')
            ylabel(ax3, 'Radius (nm)')
            out = [];

            clipboard('copy', sprintf('%.1f',median(r)))
            
            axCur = obj.initaxis('r vs curvature');
            par1 = r;
            par2 = curvature;
            ax = axCur;
            if p.comp
                hold(ax, 'on')
                for k = 1:length(grpLabel)
                    if 0
                        plot(ax, par1(fileNumber == k), par2(fileNumber == k), ' o');
                    else
                        plotSElink(ax, par1(fileNumber == k), par2(fileNumber == k),ID,se, ' o');
                    end
                    finalLabel{k} = grpLabel{k};
                end
                hold(ax, 'off')
                legend(ax,finalLabel(~cellfun(@isempty, finalLabel)))
            else
                if 0
                    plot(ax, par1, par2, ' o');
                else
                    plotSElink(ax, par1, par2,ID,se, ' o');
                end
            end
            xlabel(ax, 'Radius (nm)')
            ylabel(ax, 'Curvature (r^-1; nm^-1)')
            
            par = curvature;
            if p.comp
                axCurComp = obj.initaxis('Curvature');
                ax = axCurComp;
                if 1
                    plot(ax, par(fileNumber == 1), par(fileNumber == 3), ' o');
                else
                    plotSElink(ax, par(fileNumber == k), par(fileNumber == k),ID,se, ' o');
                end
                finalLabel = grpLabel;
%                 legend(ax,finalLabel(~cellfun(@isempty, finalLabel)))
                title(ax, 'Curvature (r^-1; nm^-1)')
                xlabel(ax, grpLabel{1})
                ylabel(ax, grpLabel{3})
            end
            
        end
        function createGrpTable(obj)
            hOld = obj.guihandles.grpTable;
            
            pos=hOld.Position;
            hNew=uitable(hOld.Parent,'Data',{},'Position',[pos(1:2)+[5 -180] 500 200]);
            obj.guihandles.grpTable = hNew;
            fileNames = getFieldAsVector(obj.locData.files.file,'name');
            [~,fileNames] = fileparts(fileNames);
            fileNames = cellstr(fileNames);
            hNew.Data = [fileNames; cell(size(fileNames))]';
            hNew.ColumnName = {'Group','Label'};
            hNew.ColumnEditable = [false true];
            hNew.ColumnWidth = {400,50};
            
            %     if isa(obj,'ROIManager.Evaluate.LocMoFitGUI')
            %         colNames={'Source', 'Rule', 'Target_fit', 'Target_usr'};
            %         % check the loaded modules and hook all the LocMoFitGUI
            %         if length(obj.locData.SE.processors.eval.guihandles.modules.Data)>1
            %             loadedModuls = obj.locData.SE.processors.eval.guihandles.modules.Data(:,2);
            %             lLocMoFitGUI = contains(loadedModuls, 'LocMoFitGUI');
            %         else
            %             lLocMoFitGUI = 0;
            %         end
            %         if sum(lLocMoFitGUI)>0
            %             loadedLocMoFitGUI = loadedModuls(lLocMoFitGUI,:);
            %         else
            %             loadedLocMoFitGUI = [];
            %         end
            %         htable.CellEditCallback = {@convertTable_callback,4};
            %         colFormat = {[{'this step'}; loadedLocMoFitGUI]',[],{'none'},[]};
            %     else
            %         colNames={'Rule', 'Target_fit', 'Target_usr'};
            %         loadedLocMoFitGUI = [];
            %         htable.CellEditCallback = {@convertTable_callback,3};
            %         colFormat = {[],{'none'},[]};
            %     end
            %     htable.ColumnName = colNames;
            %     htable.ColumnFormat = colFormat;
            %     htable.ColumnEditable = true;
            %     htable.CellSelectionCallback = {@convertTableSelection_callback,obj};
            %     htable.RowName = [];
            delete(hOld);
            
            
            
            %             h = uitab(hOld, obj);
            %             colFormat = {[],{'post_z','post_scale'},[]};
            %             h.ColumnFormat = colFormat;
            %             h.Position(3:4)=[300 150];
            %             obj.guihandles.convertTable = h;
            %             h.Data = {'find.scalingFactor','post_scale','';...
            %                 'find.scalingFactor*(150*sin(deg2rad(find.binCloseAng)))','post_z',''};
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function pard=guidef(obj)

pard.comp.object=struct('String','Comparision','Style','checkbox','Callback',{{@comp_callback,obj}});
pard.comp.position=[2,1];
pard.comp.Width=1.5;

pard.t_by.object=struct('String','By','Style','text');
pard.t_by.position=[3,1];
pard.t_by.Width=0.5;

pard.by.object=struct('String', 'file|LocMoFit','Style','popupmenu', 'Value',1);
pard.by.position=[3,1.5];
pard.by.Width=1;

pard.grpTable.object=struct('String','grpTable','Style','text');
pard.grpTable.position=[4,1];
pard.grpTable.Width=1;

% pard.comp.object=struct('String','Comparision','Style','checkbox');
% pard.comp.position=[2,1];
% pard.comp.Width=1.5;
% pard.comp.Visible='off';

%


pard.plugininfo.type='ROI_Analyze';

end

function comp_callback(a,b,obj)
if a.Value
    setVisible(obj, {'t_by' 'by' 'grpTable'}, 'on')
else
    setVisible(obj, {'t_by' 'by' 'grpTable'}, 'off')
end
end

function setVisible(obj,guiHandleName,visible)
for k = 1:length(guiHandleName)
    obj.guihandles.(guiHandleName{k}).Visible = visible;
end
end
