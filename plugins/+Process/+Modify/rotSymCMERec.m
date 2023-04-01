classdef rotSymCMERec<interfaces.DialogProcessor&interfaces.SEProcessor
% This plugin rotates the localization about the z-axis

    methods
        function obj=rotSymCMERec(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'se_viewer'};
            obj.showresults=true;
        end
        
        function out=run(obj,p)
            distFromOrigin = 2000;
            locData = obj.locData;
            pos = [];
            pos(2) = 0;
            roiSize = obj.getPar('se_siteroi');
            roiWidth = roiSize-p.spatialTrimXY(1)*2;
            sr_layerson = find(locData.getPar('sr_layerson'));
            for l = 1:length(sr_layerson)
                p_render(l) = obj.getLayerParameters(l,renderSMAP);
                p_render(l).sr_pos = [distFromOrigin+((p.binNumber+1)./2)*roiWidth 0 0];
                p_render(l).sr_size = [roiWidth*p.binNumber roiWidth]./2;
                p_render(l).sr_pixrec = 1;
            
                p_render_sideView(l) = p_render(l);
                p_render_sideView(l).sr_pos = [distFromOrigin+((p.binNumber+1)./2)*roiWidth 65 0];
                p_render_sideView(l).sr_size = [roiWidth*p.binNumber roiWidth*4/5]./2;
            end

            allAng = 0:p.rotAng:359;
            for m = 1:length(allAng)
                for k = 1:p.binNumber
                    pos(1) = k*p.distBetweenBins+distFromOrigin;
                    [locs,indLocs] = locData.getloc({'xnm','ynm', 'znm','locprecnm','layer','channel'},...
                        'grouping', 'grouped',...
                        'layer',find(locData.getPar('sr_layerson')),...
                        'position', [pos roiWidth./2],...
                        'removeFilter',{'filenumber','rank_masterAvg'}); % per ROI info.
                    locs.xnm = locs.xnm-pos(1);
                    
                    [locs.xnm,locs.ynm] = rotcoord(locs.xnm,locs.ynm,allAng(m));
                    locs.xnm = locs.xnm+rand(1)*p.jitterXY;
                    locs.ynm = locs.ynm+rand(1)*p.jitterXY;
                    
                    locs.xnm = locs.xnm+k*roiWidth+distFromOrigin;
                    if k == 1
                        newlocs = locs;
                    else
                        newlocs = mergeStruct(locs,newlocs);
                    end
                end

                for l = 1:length(sr_layerson)
                    topView_ = renderSMAP(newlocs, p_render(l), l);
                    if l == 1
                        topView = zeros(size(topView_));
                    end
                    topView_ = drawerSMAP(topView_, p_render(l));
                    topView = topView + topView_.image;
                end
                
                % get cross section
                lSection = newlocs.ynm > -30 & newlocs.ynm < 30;
                newlocs = subsetStruct(newlocs,lSection);
                
                % swap the two axes
                tempZ = newlocs.znm;
                newlocs.znm = newlocs.ynm;
                newlocs.ynm = tempZ;
                
                for l = 1:length(sr_layerson)
                    sideView_ = renderSMAP(newlocs, p_render_sideView(l), l);
                    if l == 1
                        sideView = zeros(size(sideView_));
                    end
                    sideView_ = drawerSMAP(sideView_, p_render_sideView(l));
                    sideView = sideView + sideView_.image;
                end
                if m == 1
                    finalTopView = topView;
                    finalSideView = sideView;
                else
                    finalTopView = finalTopView + topView;
                    finalSideView = finalSideView + sideView;
                end
            end
            
            fig = figure('name','RotSym');
            pan = panel(fig);
            pan.pack('v',[1 4/5]./sum([1 4/5]))
            pan(1).select()
            imagesc(finalTopView./m)
            pan(2).select()
            imagesc(finalSideView./m)
            axis(pan.de.axis, 'image')
            set(pan.de.axis, 'box','off', 'xtick',[],'ytick',[])
            pan.margin = 0;
            pan.de.margin = 0;
            addScalebar(pan(2).axis,'bottom-right',[20 20], 100)
            
            fig.Units = 'pixels';
            fig.Position(3:4) = [roiWidth*p.binNumber roiWidth*(4/5+1)]./4.01;
            exportgraphics(fig, p.filePath,'ContentType','vector', 'Resolution', 600)
            out = [];
        end
              
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function pard=guidef(obj)

pard.t_binNumber.object=struct('Style','text','String','Bin number');
pard.t_binNumber.position=[1,1];
pard.t_binNumber.Width=1;

pard.binNumber.object=struct('Style','edit','String', '10');
pard.binNumber.position=[1,2];
pard.binNumber.Width=0.5;

pard.t_distBetweenBins.object=struct('Style','text','String','Spatial bin gap');
pard.t_distBetweenBins.position=[2,1];
pard.t_distBetweenBins.Width=1;

pard.distBetweenBins.object=struct('Style','edit','String', '350');
pard.distBetweenBins.position=[2,2];
pard.distBetweenBins.Width=0.5;

pard.t_spatialTrimXY.object=struct('Style','text','String','Bin crop [X Y]');
pard.t_spatialTrimXY.position=[3,1];
pard.t_spatialTrimXY.Width=1;

pard.spatialTrimXY.object=struct('Style','edit','String', '75 75');
pard.spatialTrimXY.position=[3,2];
pard.spatialTrimXY.Width=0.5;

pard.t_masterAvgR.object=struct('Style','text','String','Radius (master avg.)');
pard.t_masterAvgR.position=[4,1];
pard.t_masterAvgR.Width=1;

pard.masterAvgR.object=struct('Style','edit','String', '150');
pard.masterAvgR.position=[4,2];
pard.masterAvgR.Width=0.5;

pard.t_rotAng.object=struct('Style','text','String','Rotaion angle');
pard.t_rotAng.position=[5,1];
pard.t_rotAng.Width=1;

pard.rotAng.object=struct('Style','edit','String','5');
pard.rotAng.position=[5,2];
pard.rotAng.Width=0.5;

pard.t_jitterXY.object=struct('Style','text','String','Jitter xy');
pard.t_jitterXY.position=[5,3];
pard.t_jitterXY.Width=1;

pard.jitterXY.object=struct('Style','edit','String','2');
pard.jitterXY.position=[5,4];
pard.jitterXY.Width=0.5;

pard.filePath.object=struct('Style','edit','String', '');
pard.filePath.position=[6,1];
pard.filePath.Width=3;

pard.saveTo.object=struct('Style','pushbutton','String', 'Save to','Callback',{{@saveTo_callback, obj}});
pard.saveTo.position=[6,4];
pard.saveTo.Width=1;


% pard.syncParameters={{'roimanager_processors','parsTable',{'Data'}}};

pard.plugininfo.description='Dynamic reconstruction of mammalian CME.';
pard.plugininfo.type='ROI_Analyze';
end

function saveTo_callback(a,b,obj)
    filter = {'*.pdf';'*.*'};
    [file,path] = uiputfile(filter);
    p = obj.getAllParameters;
    p.filePath = [path file];
    obj.setGuiParameters(p);
end