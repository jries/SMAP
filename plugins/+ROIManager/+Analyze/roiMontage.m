classdef roiMontage<interfaces.DialogProcessor&interfaces.SEProcessor
    methods
        function obj=roiMontage(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'se_viewer'};
            obj.showresults=true;
        end
        
        function out=run(obj,p)
            se = obj.locData.SE;
            numOfRoiToPlot = p.roiOrder(2)-p.roiOrder(1)+1;
            roiSize = se.P.par.se_sitefov.content;
            
            roiToPlot = cell(numOfRoiToPlot,1);
            
            %% extract individual ROIs, add their IDs
            if p.showLabel
                fig = figure('visible','off');
                a = axes(fig);
                set(fig, 'Position', [0, 0, roiSize, roiSize])
                set(a, 'Position', [0, 0, roiSize, roiSize])
                a.XLim = [0 500];
                a.YLim = [0 500];
            end
            if p.takeAll
                p.roiOrder = [1 obj.SE.numberOfSites];
            end
            for k = p.roiOrder(1):p.roiOrder(2)
                if p.onlyUsed
                    useThisRoi = se.sites(k).annotation.use;
                else
                    useThisRoi = 1;
                end
                
                if useThisRoi
                    roiToPlot{(k-p.roiOrder(1)+1),1} = se.sites(k).image.image;
                    if p.showLabel
                        % add ROIs' ID labels
                        text(a, .05,.9,num2str(se.sites(k).ID),'FontSize',38,'FontWeight','bold')
                        F = getframe(a,[0, 0, roiSize, roiSize]);
                        cla(fig)
                        F = F.cdata==0;
                        roiToPlot{(k-p.roiOrder(1)+1),1}(F==1) = 255;
                    end
                end
            end
            if p.showLabel
                close(fig)
            end
            roiToPlot = roiToPlot(~cellfun('isempty',roiToPlot));
            nrow = ceil(length(roiToPlot)/p.ncol);
            roiToPlot = cellfun(@(x)imcrop(x,[p.crop(1)/2 p.crop(2)/2 roiSize-p.crop(1) roiSize-p.crop(2)]), roiToPlot, 'UniformOutput', false);
            saveTo = [p.folder '\' p.fileName];
            img = montage(roiToPlot, 'Size', [nrow p.ncol], 'ThumbnailSize', [], 'BackgroundColor', 'white', 'BorderSize', p.pad);
            imwrite(img.CData, saveTo)
            out=[];
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end




function pard=guidef(obj)

pard.t1.object=struct('String','ROI number ([from] [to])','Style','text');
pard.t1.position=[1,1];
pard.t1.Width=1;

pard.roiOrder.object=struct('String','1 100','Style','edit', 'Enable','off');
pard.roiOrder.position=[1,2];
pard.roiOrder.Width=1;

pard.takeAll.object=struct('String','All','Style','checkbox', 'Value', 1, 'Callback', {{@setvisibility_callback,obj}});
pard.takeAll.position=[1,3];
pard.takeAll.Width=0.5;

pard.onlyUsed.object=struct('String','Only used','Style','checkbox', 'Value', 1);
pard.onlyUsed.position=[1,3.5];
pard.onlyUsed.Width=0.8;

pard.showLabel.object=struct('String','Show ID','Style','checkbox', 'Value', 1);
pard.showLabel.position=[1,4.1];
pard.showLabel.Width=0.8;

pard.t2.object=struct('String','# of column','Style','text');
pard.t2.position=[2,1];
pard.t2.Width=1;

pard.ncol.object=struct('String','5','Style','edit');
pard.ncol.position=[2,2];
pard.ncol.Width=1;

pard.t3.object=struct('String','Cropping ([x] [y])','Style','text');
pard.t3.position=[3,1];
pard.t3.Width=1;

pard.crop.object=struct('String','20 20','Style','edit');
pard.crop.position=[3,2];
pard.crop.Width=1;

pard.t4.object=struct('String','Padding ([x] [y])','Style','text');
pard.t4.position=[3,3];
pard.t4.Width=1;

pard.pad.object=struct('String','2 2','Style','edit');
pard.pad.position=[3,4];
pard.pad.Width=1;

pard.t5.object=struct('String','Save to:','Style','text');
pard.t5.position=[4,1];
pard.t5.Width=1;

pard.folder.object=struct('String','.','Style','edit');
pard.folder.position=[4,2];
pard.folder.Width=2;

pard.t6.object=struct('String','\','Style','text');
pard.t6.position=[4,4];
pard.t6.Width=0.1;

pard.fileName.object=struct('String','motage.png','Style','edit');
pard.fileName.position=[4,4.1];
pard.fileName.Width=0.9;

pard.selectFolder.object=struct('String','Select folder...','Style','pushbutton', 'Callback', {{@selectFolder_callback,obj}});
pard.selectFolder.position=[5,2];
pard.selectFolder.Width=1;





pard.plugininfo.type='ROI_Analyze';

end

function setvisibility_callback(a,b,obj)
    if obj.getSingleGuiParameter('takeAll')
        obj.guihandles.roiOrder.Enable = 'off';
    else
        obj.guihandles.roiOrder.Enable = 'on';
    end
end

function selectFolder_callback(a,b,obj)
    f = uigetdir(obj.getSingleGuiParameter('folder'));
    if f~=0
        obj.setGuiParameters(struct('folder',f));
    end
end