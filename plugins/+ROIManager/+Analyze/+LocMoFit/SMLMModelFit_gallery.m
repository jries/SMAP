classdef SMLMModelFit_gallery<interfaces.DialogProcessor&interfaces.SEProcessor
%     makes a montage of many ROIs
    % Todo: this plugin is still CME3D specific. Flexibility has to be
    % improved.
    methods
        function obj=SMLMModelFit_gallery(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'se_viewer'};
            obj.showresults=true;
        end
        
        function out=run(obj,p)
            % ask user to specify the file name
            [file,path] = uiputfile('*.png', 'Save as', '');
            file = strsplit(file,'.');
            
            % basic info.
            se = obj.locData.SE;
            sites = se.sites;
            roiSize = se.P.par.se_siteroi.content;
            
            % hack an evaluate plug-in in order to use the obj.getLocs(...)
            fdcal=figure(233);
            dcal=plugin('ROIManager','Evaluate','generalStatistics',fdcal,obj.P);
            dcal.attachLocData(obj.SE.locData);
            dcal.makeGui;
            
            % [to-do] allow selection
            fitterGUI_name = p.fitter.selection;
            eval = obj.locData.SE.processors.eval.guihandles.modules.Data(:,2);
            idxFitterGUI = strcmp(fitterGUI_name,eval);
            fitter = se.processors.eval.processors{idxFitterGUI}.fitter;
            fig = figure(234);
            fig.Position(3:4) = fig.Position(3:4).*1.2;
            pos = fig.Position;
            ax = axes(fig);
            
            page_numOfSites = 20;
            numOfView = 6;
            
            sub = cell(page_numOfSites*numOfView,1);
            subSites = se.sites(p.sites2plot);
            for k = 1:length(subSites)
                % Borrow the evaluate plug-in to use getLocs(obj,...)
                dcal.site=subSites(k);
                dcal.site.image = se.plotsite(subSites(k));
                % Check all the possible layers (up to 6)
                for l = 1:6
                    layercheck = obj.getPar(['layer' num2str(l) '_layercheck']);
                    if layercheck
                        [locsSiteOne,indlocOne] = dcal.getLocs({'xnmrot','ynmrot','znm','locprecnm', 'locprecznm'},'size',roiSize','grouping', 'grouped','layer',l); % per ROI info.
                        if l == 1
                            locsSite = locsSiteOne;
                            indloc = indlocOne;
                            locsSite.layer = ones(size(locsSiteOne.xnm));
                            locsSite.xnm = locsSite.xnmrot;
                            locsSite.ynm = locsSite.ynmrot;
                        else
                            locsSite.layer = [locsSite.layer; ones(size(locsSiteOne.xnm))*l];
                            locsSite.xnm = [locsSite.xnm; locsSiteOne.xnmrot];
                            locsSite.ynm = [locsSite.ynm; locsSiteOne.ynmrot];
                            locsSite.znm = [locsSite.znm; locsSiteOne.znm];
                            locsSite.locprecnm = [locsSite.locprecnm; locsSiteOne.locprecnm];
                            locsSite.locprecznm = [locsSite.locprecznm; locsSiteOne.locprecznm];
                        end
                    end
                end
                fitter.allParsArg = subSites(k).evaluation.(fitterGUI_name).allParsArg;
                fitter.setParArg('m1.lPar.variation', 'value',0);
                
                % [to-do] here need to generized so that the models are not
                % limited to the first one.
                fitter.model{1}.sigma = 5;
                fitter.model{1}.fixSigma = true;
                [~,modViz] = fitter.plot(locsSite,'plotType','point', 'doNotPlot', true); % get point type visualization
                lPars = fitter.exportPars(1,'lPar');
                locsViz = fitter.locsHandler(locsSite, lPars,1);
                
%                 'XYZ' is for mode
                % rotate the view to show the open
                [x,y,z] = rotcoord3(0,0,-1, deg2rad(lPars.xrot), deg2rad(lPars.yrot), deg2rad(lPars.zrot), 'XYZ');
                [aziOri,eleOri,~] = cart2sph(0,0,-1);
                [azi,ele,~] = cart2sph(x,y,z);
                
                lParsAzi = lPars;
                lParsAzi.xrot = 0;
                lParsAzi.yrot = 0;
                lParsAzi.zrot = 0;
                locsVizAzi = fitter.locsHandler(locsSite, lParsAzi,1);
                [locsVizAzi.xnm, locsVizAzi.ynm] = rotcoord2(locsVizAzi.xnm, locsVizAzi.ynm, azi);
                
                fnMod = fieldnames(modViz{1});
                for l = 1:length(fnMod)
                    modVizAzi{1}.(fnMod{l}) = modViz{1}.(fnMod{l});
                end
                [modVizAzi{1}.x, modVizAzi{1}.z] = rotcoord2(modVizAzi{1}.x, modVizAzi{1}.z, -(ele-eleOri));
                
                for rot = 1:numOfView
                    if rot == 1
                        fitter.rotCoordNMkImg(ax, modViz, locsViz, [0 0], 2, 'Data', 500, {'red hot', 'cyan cold'})
                    elseif rot == 6
                        fitter.rotCoordNMkImg(ax, modVizAzi, locsVizAzi, [0 -90], 2, 'Data', 30, {'red hot', 'cyan cold'})
                    else
                        fitter.rotCoordNMkImg(ax, modViz, locsViz, [45*(rot-2) -90], 2, 'Data', 30, {'red hot', 'cyan cold'})
                    end
                    ax.Children(1).Color = [0.7 0.7 0.7];
                    set(ax,'YDir','normal')
                    if rot == 1
                        text(ax, page_numOfSites, 30, num2str(subSites(k).ID),'FontSize',50, 'Color','w')
                    end
                    currentFrame = getframe(fig);
                    mon_order = rem(k,page_numOfSites);
                    if mon_order == 0
                        mon_order = page_numOfSites;
                    end
                    sub{(mon_order-1)*numOfView+rot} = currentFrame.cdata(40:440,144:544,:);
                end
                if rem(k,page_numOfSites)==0 || k==length(subSites)
                    imgm = montage(sub, 'Size', [page_numOfSites numOfView], 'ThumbnailSize', [], 'BackgroundColor', 'white', 'BorderSize', 2);
                    imwrite(imgm.CData, [path file{1} '_' num2str(ceil(k/page_numOfSites)) '.' file{2}])
                    close(fig)
                    fig = figure(234);
                    ax = axes(fig);
                    fig.Position = pos;
                    sub = [];
                end
            end
            close(fig)
            out = [];
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end




function pard=guidef(obj)

pard.t1.object=struct('String','Which fitter','Style','text');
pard.t1.position=[1,1];
pard.t1.Width=1;

eval = obj.locData.SE.processors.eval.guihandles.modules.Data(:,2);
lFitter = startsWith(eval,'SMLMModelFitGUI');
options = eval(lFitter);

pard.fitter.object=struct('String',{options},'value',1,'Style','popupmenu');
pard.fitter.position=[1,2];
pard.fitter.Width=1;

pard.t2.object=struct('String','Sites to plot','Style','text');
pard.t2.position=[2,1];
pard.t2.Width=1;

pard.sites2plot.object=struct('String', 'Sites Index','Style','edit');
pard.sites2plot.position=[2,2];
pard.sites2plot.Width=1;

pard.plugininfo.description='makes a montage of many ROIs';
pard.plugininfo.type='ROI_Analyze';
end
