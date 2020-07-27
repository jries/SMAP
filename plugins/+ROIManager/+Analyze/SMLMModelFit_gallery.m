classdef SMLMModelFit_gallery<interfaces.DialogProcessor&interfaces.SEProcessor
%     makes a montage of many ROIs
    methods
        function obj=SMLMModelFit_gallery(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'se_viewer'};
            obj.showresults=true;
        end
        
        function out=run(obj,p)
            [file,path] = uiputfile('*.png', 'Save as', '');
            file = strsplit(file,'.');
            se = obj.locData.SE;
            sites = se.sites;
            roiSize = se.P.par.se_siteroi.content;
            
            % hack an evaluate plug-in in order to use the obj.getLocs(...)
            fdcal=figure(233);
            dcal=plugin('ROIManager','Evaluate','generalStatistics',fdcal,obj.P);
            dcal.attachLocData(obj.SE.locData);
            dcal.makeGui;
            
            % to-do: allow selection
            fitterGUI_name = p.fitter.selection;
            eval = obj.locData.SE.processors.eval.guihandles.modules.Data(:,2);
            idxFitterGUI = find(strcmp(fitterGUI_name,eval));
            fitter = se.processors.eval.processors{idxFitterGUI}.fitter;
            fig = figure(234);
            fig.Position(3:4) = fig.Position(3:4).*1.2;
            pos = fig.Position;
            ax = axes(fig);
            sub = cell(20*5,1);
            % Borrow the evaluate plug-in to use getLocs(obj,...)
            for k = 1:se.numberOfSites
                dcal.site=sites(k);
                dcal.site.image = se.plotsite(sites(k));
                [locsSite,indloc] = dcal.getLocs({'xnmrot','ynmrot','znm','locprecnm', 'locprecznm'},'size',roiSize','grouping', 'grouped','layer',1); % per ROI info.
                locsSite.layer = ones(size(locsSite.xnm));
                locsSite.xnm = locsSite.xnmrot;
                locsSite.ynm = locsSite.ynmrot;
                fitter.allParsArg = sites(k).evaluation.SMLMModelFitGUI.allParsArg;
                fitter.setParArg('m1.lPar.variation', 'value',0);
                [~,modViz] = fitter.plot(locsSite,'plotType','point', 'doNotPlot', true); % get point type visualization
                locsViz = fitter.locsHandler(locsSite,fitter.exportPars(1,'lPar'),1);
                
                for rot = 1:5
                    if rot == 1
                        fitter.rotCoordNMkImg(ax, modViz, locsViz, [0 0], 2, 'Data', 500, {'red hot'})
                    else
                        fitter.rotCoordNMkImg(ax, modViz, locsViz, [45*(rot-2) -90], 2, 'Data', 30, {'red hot'})
                    end
                    ax.Children(1).Color = [0.7 0.7 0.7];
                    set(ax,'YDir','normal')
                    if rot == 1
                        text(ax, 20, 30, num2str(sites(k).ID),'FontSize',50, 'Color','w')
                    end
                    currentFrame = getframe(fig);
                    mon_order = rem(k,20);
                    if mon_order == 0
                        mon_order = 20;
                    end
                    sub{(mon_order-1)*5+rot} = currentFrame.cdata(40:440,144:544,:);
                end
                if rem(k,20)==0 || k==se.numberOfSites
                    imgm = montage(sub, 'Size', [20 5], 'ThumbnailSize', [], 'BackgroundColor', 'white', 'BorderSize', 2);
                    imwrite(imgm.CData, [path file{1} '_' ceil(num2str(k/20)) '.' file{2}])
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

pard.t1.object=struct('String','Which fitter:','Style','text');
pard.t1.position=[1,1];
pard.t1.Width=1;

eval = obj.locData.SE.processors.eval.guihandles.modules.Data(:,2);
lFitter = startsWith(eval,'SMLMModelFitGUI');
options = eval(lFitter);

pard.fitter.object=struct('String',options,'Value',1,'Style','popupmenu');
pard.fitter.position=[1,2];
pard.fitter.Width=1;

pard.plugininfo.description='makes a montage of many ROIs';
pard.plugininfo.type='ROI_Analyze';
end
