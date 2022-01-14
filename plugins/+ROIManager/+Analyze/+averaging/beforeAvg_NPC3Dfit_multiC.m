classdef beforeAvg_NPC3Dfit_multiC<interfaces.DialogProcessor&interfaces.SEProcessor
    methods
        function obj=beforeAvg_NPC3Dfit_multiC(varargin)        
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
            fitInfo = getFieldAsVector(usedSites, 'evaluation.SMLMModelFitGUI_2.fitInfo');
            lFailed = cellfun(@(x)strcmp(x.guiInfo,'Fit or plot failed.'), fitInfo);
            evalList = se.processors.eval.guihandles.modules.Data(:,2);
            indProcessor = find(strcmp('SMLMModelFitGUI',evalList));
            relativePosLastStep = 1;
            
            ID = getFieldAsVector(usedSites, 'ID');
            numOfLocsL1 = getFieldAsVectorInd(usedSites, 'evaluation.SMLMModelFitGUI_2.fitInfo.numOfLocsPerLayer', 1);
            % [~,idxZ] = g.locData.SE.processors.eval.processors{indProcessor-2}.fitter.wherePar('pars.m1.lPar.z');
            % zS1 = getFieldAsVectorInd(usedSites, 'evaluation.SMLMModelFitGUI.allParsArg.value',idxZ);
            % 
            
            [~,idxRingDist] = se.processors.eval.processors{indProcessor}.fitter.wherePar('pars.m2.lPar.z');
            ringDistS1 = getFieldAsVectorInd(usedSites, 'evaluation.SMLMModelFitGUI.allParsArg.value',idxRingDist);

%             [cutoffOneRing, bin_edges] = getOneRingCutoff(ringDistS1);
            
            [~,idxZ] = se.processors.eval.processors{indProcessor+relativePosLastStep}.fitter.wherePar('pars.m1.lPar.z');
            z = getFieldAsVectorInd(usedSites, 'evaluation.SMLMModelFitGUI_2.allParsArg.value',idxZ);
            
            lOneRing = ringDistS1<=cutoffOneRing;
            
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
            
                        
            %% Sort sites
            %%
            g = obj.getPar('mainGui');
            
            sortROIs = g.children.guiSites.children.Helper.children.SortROIs;
            guiPar = sortROIs.getGuiParameters;
            guiPar.direction1.Value = 2;
            guiPar.sortedit1 = 'annotation.use';
            guiPar.sortprop1.Value = 3;

            guiPar.direction2.Value = 1;
            guiPar.sortedit2 = 'annotation.list3.value';
            guiPar.sortprop2.Value = 2;
            
            % disable all other sorts
            guiPar.sortedit3 = '';
            guiPar.sortprop3.Value = 1;
            guiPar.sortedit4 = '';
            guiPar.sortprop4.Value = 1;

            sortROIs.setGuiParameters(guiPar);
            sortROIs.run(sortROIs.getAllParameters);
            sites = se.sites;
            
            %% Clear alignment
            if p.reAverage
                fn = fieldnames(obj.locData.loc);
                lRm = contains(fn, '_SMLMModelFitGUI_');
                if ~isempty(lRm)
                    obj.locData.loc = rmfield(obj.locData.loc, fn(lRm));
                    obj.locData.grouploc = rmfield(obj.locData.grouploc, fn(lRm));
                end
                
                %% Redraw sites
                seGUI = g.children.guiSites.SEpreview;
                list3 = getFieldAsVector(sites, 'annotation.list3.value');
                lUse = getFieldAsVector(sites, 'annotation.use');
                lastSite = find(list3==1&lUse,1,'last');
                sitelist = seGUI.getSingleGuiParameter('sitelist');
                sitelist.Value = 1:lastSite;
                seGUI.setGuiParameters(struct('sitelist', sitelist));
                
                % set 
                fitStep1 = se.processors.eval.children.SMLMModelFitGUI;
                guiPar = fitStep1.getGuiParameters;
                guiPar.noFit = 1;
                guiPar.useAlignment = 0;
                fitStep1.setGuiParameters(guiPar);
                
                fitStep2 = se.processors.eval.children.SMLMModelFitGUI_2;
                guiPar = fitStep2.getGuiParameters;
                guiPar.noFit = 1;
                guiPar.useAlignment = 1;
                fitStep2.setGuiParameters(guiPar);
                
                seGUI.redrawall;
                se.locData.regroup;
                
                %% Create avg.
                % define main coord
                module=plugin('Process','Modify','DefineMainCoordinates');
                p.Vrim=100;

                module.handle=figure('MenuBar','none','Toolbar','none','Name','DefineMainCoordinates');
                module.attachPar(obj.P);
                module.attachLocData(obj.locData);

                p.Xrim=10;
                module.setGuiAppearence(p)
                module.makeGui;

                guiPar = module.getGuiParameters;
                guiPar.changefield = 1; guiPar.changex = 1; guiPar.changey = 1; guiPar.changez = 1;
                allOptions = guiPar.xfield.String;
                idx = find(strcmp(allOptions, 'xnmaligned_SMLMModelFitGUI_2'));
                guiPar.xfield.Value = idx;
                idx = find(strcmp(allOptions, 'ynmaligned_SMLMModelFitGUI_2'));
                guiPar.yfield.Value = idx;
                idx = find(strcmp(allOptions, 'znmaligned_SMLMModelFitGUI_2'));
                guiPar.zfield.Value = idx;
                module.setGuiParameters(guiPar);
                guiPar = module.getGuiParameters;
                module.run(guiPar);
            
            end
            
            
            
%             L1 = g.children.guiRender.children.Layer1;
%             guiPar = L1.getGuiParameters;
%             guiPar.colorfield_min = -150;
%             guiPar.colorfield_max = 150;
%             guiPar.imax_min = -2.4;
%             guiPar.lut.Value = 4;
%             guiPar.render_colormode.Value = 2;
%             L1.setGuiParameters(guiPar)
% 
%             guiFormat = g.getPar('guiFormat');
%             guiFormat.guihandles.pixrec.String = '1';
%             g.setPar('sr_pixrec',3)
%             g.setPar('sr_colorbarthickness',40)
%             g.setPar('sr_plotscalebar',0)
%             g.setPar('sr_pos',[34391 47036])
% 
%             sr_fig = g.getPar('sr_figurehandle');
%             sr_fig.Position(3:4) = [760 730];
% 
%             notify(g.P,'sr_render')

            
            
            %% Clear alignment
            if p.saveFile
                % rm locs
                L1 = g.children.guiRender.children.Layer1;
                guiPar = L1.getGuiParameters;
                guiPar.znm_min = -130;
                guiPar.znm_max = 130;
                guiPar.imax_min = -3.5;
                guiPar.lut.Value = 10;
                guiPar.lutinv = 0;
                guiPar.layercheck = 1;
                L1.setGuiParameters(guiPar)
                L1.updateLayerField

                L2 = g.children.guiRender.children.Layer2;
                guiPar = L2.getGuiParameters;
                guiPar.znm_min = -130;
                guiPar.znm_max = 130;
                guiPar.imax_min = -3.5;
                guiPar.lut.Value = 9;
                guiPar.lutinv = 0;
                guiPar.layercheck = 1;
                L2.setGuiParameters(guiPar)
                L2.updateLayerField

    %             guiFormat = g.getPar('guiFormat');
    %             guiFormat.guihandles.pixrec.String = '1';
                g.setPar('sr_pixrec',0.7)
                g.setPar('sr_pos',[500 500])
                notify(g.P,'sr_render')
                
                module=plugin('Process','Modify','RemoveLocs');
                p.Vrim=100;

                module.handle=figure('MenuBar','none','Toolbar','none','Name','RemoveLocs');
                module.attachPar(obj.P);
                module.attachLocData(obj.locData);

                p.Xrim=10;
                module.setGuiAppearence(p)
                module.makeGui;

                guiPar = module.getGuiParameters;
                guiPar.roic.Value = 3;
                guiPar.allfiles = 1;
                module.setGuiParameters(guiPar);
                module.run(module.getGuiParameters);

                % impose 8-fold sym
                module=plugin('Process','Modify','ImposeSymmetry');
                p.Vrim=100;

                module.handle=figure('MenuBar','none','Toolbar','none','Name','ImposeSymmetry');
                module.attachPar(obj.P);
                module.attachLocData(obj.locData);

                p.Xrim=10;
                module.setGuiAppearence(p)
                module.makeGui;
                module.run(module.getGuiParameters);

                %% Save file
                oldName = g.locData.files.file.name;

                today = char(datetime('today','Format','yyMMdd'));
                originalDate = regexp(oldName, '\_c(\d{6})', 'tokens');
                originalDate = char(originalDate{:}{:});
                lRm = originalDate == today;

                newName = regexprep(oldName, '\_d(\d){0,6}\_', '');
                newName = regexprep(newName, '\_c(\d){6}', ['_d' originalDate(~lRm)]);
                newName = [newName(1:end-8) ['_c' today '_avg8FImp'] newName(end-7:end)];
                g.locData.savelocs(newName);
            end
            %
            %%
            %% Ring separation
            axRS=obj.initaxis('Ring separation vs z');
            f = fit(z(lGood),ringDistS1(lGood),'poly1','Robust','LAR');
            h1 = plot(axRS, z(lGood),ringDistS1(lGood), ' ob');
            hold(axRS,'on')
            h2 = plot(axRS, z(~lGood),ringDistS1(~lGood), ' ok');
            h3 = plot(f, 'c');
            hold(axRS,'off')
            xlabel(axRS, 'z position (nm)')
            ylabel(axRS, 'Ring separation (nm)')
            legend([h1 h2 h3],{'Included','Excluded','Linear fit'})           

            out = [];
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function pard=guidef(obj)

pard.reAverage.object=struct('String','Re-run average','Style','checkbox','Value',0);
pard.reAverage.position=[2,1];
pard.reAverage.Width=1.5;

pard.saveFile.object=struct('String','Save','Style','checkbox','Value',0);
pard.saveFile.position=[3,1];
pard.saveFile.Width=1.5;

pard.plugininfo.type='ROI_Analyze';

end