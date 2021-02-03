classdef siteGallery_mCME<interfaces.DialogProcessor&interfaces.SEProcessor
    methods
        function obj=siteGallery_mCME(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'se_viewer', 'se_siteroi', 'se_sitefov'};
            obj.showresults=true;
        end
        
        function out=run(obj,p)
            global pan
            % basic info.
            se = obj.locData.SE;
            roiSize = se.P.par.se_siteroi.content;
            sites = se.sites;
            sites2plot = p.sites;
            
            % hack an evaluate plug-in in order to use the obj.getLocs(...)
            fdcal=figure(233);
            dcal=plugin('ROIManager','Evaluate','generalStatistics',fdcal,obj.P);
            dcal.attachLocData(obj.SE.locData);
            dcal.makeGui;
            
            % [to-do] allow selection
            fitterGUI_name = 'SMLMModelFitGUI';
            eval = obj.locData.SE.processors.eval.guihandles.modules.Data(:,2);
            idxFitterGUI = strcmp(fitterGUI_name,eval);
            fitter = copy(se.processors.eval.processors{idxFitterGUI}.fitter);
            
            obj.setPar('fitter',fitter)
            
            numOfView = 2;
            
            siteID = getFieldAsVector(se.sites, 'ID');
            
            [~,siteOrder] = ismember(sites2plot,siteID);
            subSites = se.sites(siteOrder);
            f=figure;
            pan = panel(f);
            pan.pack('v', {2/3 1/3});
            pan(1).pack(2, length(subSites));
            pan(2).pack(1, length(subSites));

            numOfPickedSites = length(subSites);
            obj.setPar('numOfPickedSites',numOfPickedSites);
            obj.setPar('numOfPickedSites',numOfPickedSites);
            
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
                fitter.model{1}.sigma = p.isoBlurr;
                fitter.model{1}.fixSigma = true;
                fitter.model{1}.pixelSize = p.pixelSize;
                fitter.refPoint_spacing = p.isoGap;
                fitter.roiSize = 500;
                [~,modViz] = fitter.plot(locsSite,'plotType','point', 'doNotPlot', true); % get point type visualization
                lPars = fitter.exportPars(1,'lPar');
                locsViz = fitter.locsHandler(locsSite, lPars,1);
                
%                 'XYZ' is for model
                % rotate the view to show the open
                
                
                for rot = 1:numOfView
                    if rot == 1
                        ax = fitter.rotCoordNMkImg(modViz, locsViz, [0 0], 2, 'Data', 500, {'red hot', 'cyan cold'});
                    else
                        ax = fitter.rotCoordNMkImg(modViz, locsViz, [45*(rot-2) -90], 2, 'Data', 30, {'red hot', 'cyan cold'});
                    end
                    set(ax,'YDir','normal')
%                     if rot == 1
%                        text(ax, page_numOfSites, 30, num2str(subSites(k).ID),'FontSize',50, 'Color','w')
%                     end
                    set(ax,'XTick',[], 'YTick', []);
                    axis(ax, 'image')
                    tempFig = ax.Parent;
                    pan(1,rot,k).select(ax);
                    close(tempFig);
                end
                ax = fitter.model{1}.patchPlot(fitter.exportPars(1,'mPar'));
                set(ax,'XLim',[0 fitter.roiSize/p.pixelSize]);
                set(ax,'XTick',[], 'YTick', [], 'ZTick', []);
                set(ax, 'XColor', 'none', 'YColor', 'none', 'ZColor', 'none');
                tempFig = ax.Parent;
                pan(2,1,k).select(ax);
                close(tempFig);
            end
            update_callback([],[],obj);
            out = [];
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function update_callback(a,b,obj)
    global pan
    
    p = obj.getAllParameters;
    fitter = obj.getPar('fitter');
    roiSize = fitter.roiSize;
    pixelSize = p.pixelSize;
    numOfPickedSites = obj.getPar('numOfPickedSites');
    
    pan.de.margin = 0;
    pan.margin = [2 2 2 2];
    
    f = pan.figure;
    f.Units = 'centimeters';
    f.Position(3:4) = [numOfPickedSites*4*(roiSize-p.crop*2)/roiSize+0.4 3*4*(roiSize-p.crop*2)/roiSize+0.4];
    p1 = pan(1).de.object;
    p2 = pan(2).de.object;
    p1Line = findobj(p1,'Type','line');
    axis(p1, 'image')
    axis(p2, 'equal')
%     axis(p2, 'square')
    set(p1, 'XLim', [p.crop/pixelSize (roiSize-p.crop)/pixelSize]);
    set(p1, 'YLim', [p.crop/pixelSize (roiSize-p.crop)/pixelSize]);
    set(p2, 'XLim', [p.crop/pixelSize (roiSize-p.crop)/pixelSize]);
    set(p1Line, 'Color', p.lineColor)
    set(p1Line, 'LineWidth', p.lineWidth)
end

function pard=guidef(obj)

rowRun = 1;
pard.t1.object=struct('String','Site IDs','Style','text');
pard.t1.position=[rowRun,1];
pard.t1.Width=1;

pard.sites.object=struct('String','','Style','edit');
pard.sites.position=[rowRun,2];
pard.sites.Width=1;

pard.t_pixelSize.object=struct('String','Pixel size','Style','text');
pard.t_pixelSize.position=[rowRun+1,1];
pard.t_pixelSize.Width=1;

pard.pixelSize.object=struct('String','2','Style','edit');
pard.pixelSize.position=[rowRun+1,2];
pard.pixelSize.Width=1;

pard.t_isoBlurr.object=struct('String','isoBlurr','Style','text');
pard.t_isoBlurr.position=[rowRun+2,1];
pard.t_isoBlurr.Width=1;

pard.isoBlurr.object=struct('String','15','Style','edit');
pard.isoBlurr.position=[rowRun+2,2];
pard.isoBlurr.Width=1;
pard.isoBlurr.Tooltip = 'The gaussian sigma for blurring. This controls the smoothness of the isosurface model rendering.';

pard.t_isoGap.object=struct('String','isoGap','Style','text');
pard.t_isoGap.position=[rowRun+3,1];
pard.t_isoGap.Width=1;

pard.isoGap.object=struct('String','0.15','Style','edit');
pard.isoGap.position=[rowRun+3,2];
pard.isoGap.Width=1;
pard.isoGap.Tooltip = 'Gap between sampled points. This controls the sampling rate of the isosurface model rendering.';

rowUpdate = 6;

pard.update.object=struct('String','update','Style','pushbutton','callback',{{@update_callback,obj}});
pard.update.position=[rowUpdate+0.5,3.5];
pard.update.Width=1;

pard.t_crop.object=struct('String','Crop','Style','text');
pard.t_crop.position=[rowUpdate,1];
pard.t_crop.Width=1;

pard.crop.object=struct('String','0','Style','edit');
pard.crop.position=[rowUpdate,2];
pard.crop.Width=1;

pard.t_lineColor.object=struct('String','Line color','Style','text');
pard.t_lineColor.position=[rowUpdate+1,1];
pard.t_lineColor.Width=1;

pard.lineColor.object=struct('String','#0247ff','Style','edit');
pard.lineColor.position=[rowUpdate+1,2];
pard.lineColor.Width=1;

pard.t_lineWidth.object=struct('String','Line width','Style','text');
pard.t_lineWidth.position=[rowUpdate+2,1];
pard.t_lineWidth.Width=1;

pard.lineWidth.object=struct('String','2.5','Style','edit');
pard.lineWidth.position=[rowUpdate+2,2];
pard.lineWidth.Width=1;



% pard.t2.object=struct('String','Which fitter','Style','text');
% pard.t2.position=[3,1];
% pard.t2.Width=1;
% 
% eval = obj.locData.SE.processors.eval.guihandles.modules.Data(:,2);
% lFitter = startsWith(eval,'SMLMModelFitGUI');
% options = eval(lFitter);
% 
% pard.fitter.object=struct('String',{options},'value',1,'Style','popupmenu');
% pard.fitter.position=[3,2];
% pard.fitter.Width=1;

pard.plugininfo.type='ROI_Analyze';

end