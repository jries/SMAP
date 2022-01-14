classdef summarize_yCME<interfaces.DialogProcessor&interfaces.SEProcessor
    methods
        function obj=summarize_yCME(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'se_viewer', 'se_siteroi', 'se_sitefov'};
            obj.showresults=true;
        end
        
        function out=run(obj,p)
            se = obj.locData.SE;
            sites = se.sites;

            % get used sites
            lUsed = getFieldAsVector(sites, 'annotation.use');
            usedSites = sites(lUsed);
            
            fn = fieldnames(se.processors.eval.children);
            idxFitter = find(strcmp('SMLMModelFitGUI',fn));
            fitter = se.processors.eval.processors{idxFitter}.fitter;
            
            % kick the wide sites out
            wideOrNot = getFieldAsVector(usedSites, 'annotation.list4.value');
            
            
%             list1 = getFieldAsVector(usedSites, 'annotation.list2.value');
%             list2 = getFieldAsVector(usedSites, 'annotation.list3.value');
            
 
%             lGood = wideOrNot==1;
%             goodSites = usedSites(lGood);
            
            [~,idxAxialLength] = fitter.getVariable('m2.b');
            axialLength = getFieldAsVectorInd(usedSites, 'evaluation.SMLMModelFitGUI.allParsArg.value',idxAxialLength);
            
            [~,idxBG_l1] = fitter.getVariable('m91.weight');
            BG_l1 = getFieldAsVectorInd(usedSites, 'evaluation.SMLMModelFitGUI.allParsArg.value',idxBG_l1);
            
            [~,idxBG_l2] = fitter.getVariable('m92.weight');
            BG_l2 = getFieldAsVectorInd(usedSites, 'evaluation.SMLMModelFitGUI.allParsArg.value',idxBG_l2);
            
            numOfLocs_layer1 = getFieldAsVectorInd(usedSites, 'evaluation.SMLMModelFitGUI.fitInfo.numOfLocsPerLayer',1);
            numOfLocs_layer2 = getFieldAsVectorInd(usedSites, 'evaluation.SMLMModelFitGUI.fitInfo.numOfLocsPerLayer',2);
            
            bgLocs_L1 = BG_l1.*numOfLocs_layer1;
            bgLocs_L2 = BG_l2.*numOfLocs_layer2;
            
            % based on the median numbers of the backgrounds across the
            % entire dataset, determine which sites have no meaninful
            % information.
            lNoLas17 = numOfLocs_layer1<median(bgLocs_L1);
            lNoAbp1 = numOfLocs_layer2<median(bgLocs_L2);
                       
            for k = 1:length(sites)
                sites(k).annotation.list3.value = 1;
            end
            
            for k = find(lNoAbp1)'
                usedSites(k).annotation.list3.value = 2;
            end

            for k = find(lNoLas17)'
                usedSites(k).annotation.list3.value = 3;
            end

            
            numOfLocs_BF_layer1 = getFieldAsVectorInd(usedSites, 'evaluation.SMLMModelFitGUI.fitInfo.numOfLocsPerLayer_BGFree',1);
            numOfLocs_BF_layer2 = getFieldAsVectorInd(usedSites, 'evaluation.SMLMModelFitGUI.fitInfo.numOfLocsPerLayer_BGFree',2);
                       
            g = obj.getPar('mainGui');
            sortROIs = g.children.guiSites.children.Helper.children.SortROIs;
            sortROIs.guihandles.direction1.Value = 2;
            sortROIs.guihandles.sortedit1.String = 'annotation.list2.value';
            sortROIs.guihandles.sortprop1.Value = 3;

            sortROIs = g.children.guiSites.children.Helper.children.SortROIs;
            sortROIs.guihandles.direction2.Value = 2;
            sortROIs.guihandles.sortedit2.String = 'annotation.list4.value';
            sortROIs.guihandles.sortprop2.Value = 3;
            
            sortROIs = g.children.guiSites.children.Helper.children.SortROIs;
            sortROIs.guihandles.direction3.Value = 2;
            sortROIs.guihandles.sortedit3.String = 'annotation.list3.value';
            sortROIs.guihandles.sortprop3.Value = 3;
            
            sortROIs = g.children.guiSites.children.Helper.children.SortROIs;
            sortROIs.guihandles.direction4.Value = 1;
            sortROIs.guihandles.sortedit4.String = 'evaluation.SMLMModelFitGUI.allParsArg.value(23)';
            sortROIs.guihandles.sortprop4.Value = 4;
            
            sortROIs.run(sortROIs.getAllParameters);
            
            sites = se.sites;
            lUsed = getFieldAsVector(sites, 'annotation.use');
            usedSites = sites(lUsed);
            
            wideOrNot = getFieldAsVector(usedSites, 'annotation.list4.value');
            withProOrNot = getFieldAsVector(usedSites, 'annotation.list3.value');
            goodOrBad = getFieldAsVector(usedSites, 'annotation.list2.value');
            lGood = wideOrNot==1&withProOrNot<=2&goodOrBad==1;
            
            siteOrder = 1:sum(lGood);
            axialLength = getFieldAsVectorInd(usedSites, 'evaluation.SMLMModelFitGUI.allParsArg.value',idxAxialLength);
            BG_l1 = getFieldAsVectorInd(usedSites, 'evaluation.SMLMModelFitGUI.allParsArg.value',idxBG_l1);
            BG_l2 = getFieldAsVectorInd(usedSites, 'evaluation.SMLMModelFitGUI.allParsArg.value',idxBG_l2);
            numOfLocs_layer1 = getFieldAsVectorInd(usedSites, 'evaluation.SMLMModelFitGUI.fitInfo.numOfLocsPerLayer',1);
            numOfLocs_layer2 = getFieldAsVectorInd(usedSites, 'evaluation.SMLMModelFitGUI.fitInfo.numOfLocsPerLayer',2);
            
            bgLocs_L1 = BG_l1.*numOfLocs_layer1;
            bgLocs_L2 = BG_l2.*numOfLocs_layer2;
            pseudotime = (siteOrder-1)./(max(siteOrder)-1);
            pseudotime = [zeros(1,find(lGood,1)-1)-1 pseudotime];
            
            indFirstSite = find(lGood,1);
            
            xx = getHistogramEdge(pseudotime,0.1);
            yy2_l1=bindata(pseudotime,numOfLocs_BF_layer1,xx,'median');
%             yy3_l1=bindata(1:length(numOfLocs_BF_layer1),numOfLocs_BF_layer1,xx,'robustmean');
            yy2_axialLength=bindata(pseudotime,axialLength,xx,'median');
            f_axialLen = fit(pseudotime(indFirstSite:end)',axialLength(indFirstSite:end),'poly1', 'Robust','LAR');
            
            yy2_l2=bindata(pseudotime,numOfLocs_BF_layer2,xx,'median');
%             yy3_l2=bindata(1:length(numOfLocs_BF_layer1),numOfLocs_BF_layer2,xx,'robustmean');
%             LLfit = [];
            f_l2 = fit(pseudotime(indFirstSite:end)',numOfLocs_BF_layer2(indFirstSite:end),'a.*(log(c*x+1)/log(b))','Lower',[0,1,0],'Upper',[Inf,Inf,Inf], 'Robust','LAR');
            
            axLayer1 = obj.initaxis('Locs_layer1');
            plot(axLayer1, pseudotime,numOfLocs_BF_layer1,' o')
            hold(axLayer1, 'on')
            plot(axLayer1, xx,yy2_l1)
            yline(axLayer1, median(numOfLocs_BF_layer1))
            hold(axLayer1, 'off')
            xlabel('Pseudotime (a.u.)')
            ylabel('Count (locs)')
            
            axLayer2 = obj.initaxis('Locs_layer2');
            plot(axLayer2, pseudotime,numOfLocs_BF_layer2,' o')
            hold(axLayer2, 'on')
            plot(axLayer2, xx,yy2_l2)
            plot(f_l2)
            hold(axLayer2, 'off')
            xlabel('Pseudotime (a.u.)')
            ylabel('Count (locs)')
            
            axBG_l1 = obj.initaxis('BG_layer1');
            plot(axBG_l1, pseudotime,bgLocs_L1,' o')
            hold(axBG_l1, 'on')
%             plot(axBG_l2, xx,yy2_l2)
            yline(axBG_l1, median(bgLocs_L1))
%             plot(f_l2)
            hold(axBG_l1, 'off')
            xlabel('Pseudotime (a.u.)')
            ylabel('Background locs')
            
            axBG_l2 = obj.initaxis('BG_layer2');
            plot(axBG_l2, pseudotime(indFirstSite:end),bgLocs_L2(indFirstSite:end),' ob')
            hold(axBG_l2, 'on')
%             plot(axBG_l2, xx,yy2_l2)
            plot(axBG_l2, pseudotime(1:indFirstSite-1),bgLocs_L2(1:indFirstSite-1),' ok')
            yline(axBG_l2, median(BG_l2(indFirstSite:end)))
%             plot(f_l2)
            hold(axBG_l2, 'off')
            xlabel('Pseudotime (a.u.)')
            ylabel('Background locs')
            
            axAxialLength = obj.initaxis('Axial length');
            plot(axAxialLength, pseudotime(indFirstSite:end),axialLength(indFirstSite:end),' ob')
            hold(axAxialLength, 'on')
            plot(axAxialLength, pseudotime(1:indFirstSite-1),axialLength(1:indFirstSite-1),' ok')
            plot(axAxialLength, xx,yy2_axialLength)
            plot(f_axialLen)
            hold(axAxialLength, 'off')
            xlabel('Pseudotime (a.u.)')
            ylabel('Abp1 axial length (nm)')
            
            locs_l1 = num2str(median(numOfLocs_BF_layer1));
            locs_l2 = [num2str(f_l2.a) '.*(log(' num2str(f_l2.c) '.*rand+1)./log(' num2str(f_l2.b) '))'];
            BG_l1 = num2str(median(BG_l1));
            BG_l2 = num2str(median(BG_l2(indFirstSite:end)));
            axialLen = [num2str(f_axialLen.p1) '*rand+' num2str(f_axialLen.p2)];
            
            clipboard('copy',['locs_l1:' locs_l1 newline...
                'locs_l2:' locs_l2 newline...
                'BG_l1: ' BG_l1 newline...
                'BG_l2:' BG_l2 newline...
                'axialLen:' axialLen])
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

pard.showIntActPlot.object=struct('String','Display interactive','Style','checkbox');
pard.showIntActPlot.position=[2,1];
pard.showIntActPlot.Width=1.5;

pard.showExampleMod.object=struct('String','Display example model','Style','checkbox');
pard.showExampleMod.position=[3,1];
pard.showExampleMod.Width=1.5;

pard.plugininfo.type='ROI_Analyze';

end