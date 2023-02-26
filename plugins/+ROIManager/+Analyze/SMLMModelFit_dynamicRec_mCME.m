classdef SMLMModelFit_dynamicRec_mCME<interfaces.DialogProcessor&interfaces.SEProcessor
% Export the results of the SMLMModelFitGUI
    properties
        linkedManager
    end
    properties (Dependent)
        variableIDs
    end
    methods
        function obj=SMLMModelFit_dynamicRec_mCME(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'se_viewer'};
            obj.showresults=true;
        end
        
        function makeGui(obj,varargin)
            makeGui@interfaces.DialogProcessor(obj); %make the main GUI from the guidef definitions
                %Settings
            obj.createConvertTable;
        end
        
        function out=run(obj,p)
            fit_manager = obj.linkedManager;
            
            fit_manager.parentObj.loadData;
            fit_manager.handles.linkedDynamicRec = obj;
            
            % filter the sites based on the different parameters
            boundCurvature = [-inf inf]; % take all sites
            if p.withoutClouds
                boundCurvature(2) = p.cloudsThreshold;
            end
            if p.noNeg
                boundCurvature(1) = 0;
            end
            % site filtering on sites
            fit_manager.addFiltering('LocMoFitGUI_2.m1.curvature', boundCurvature)
            fit_manager.dynamicRec;
            obj.locData.regroup;
            obj.locData.filter;
            out = [];
        end
        
        function out=createConvertTable(obj,p)
            %
            hOld = obj.guihandles.convertTable;
            h = createConvertTable(hOld, obj);
            colFormat = {[],{'post_z','post_scale'},[]};
            h.ColumnFormat = colFormat;
            h.Position(3:4)=[300 150];
            obj.guihandles.convertTable = h;
            h.Data = {'find.scalingFactor','post_scale','';...
                'find.scalingFactor*(150*sin(deg2rad(find.binCloseAng)))','post_z',''};
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

pard.t_lastBin_noRescale.object=struct('Style','text','String','No rescale up to bin');
pard.t_lastBin_noRescale.position=[5,1];
pard.t_lastBin_noRescale.Width=1;

pard.lastBin_noRescale.object=struct('Style','edit','String','1');
pard.lastBin_noRescale.position=[5,2];
pard.lastBin_noRescale.Width=0.5;

pard.noNeg.object=struct('Style','checkbox','Value',0,'String','Only curvature>0');
pard.noNeg.position=[1,3];
pard.noNeg.Width=1;

pard.withoutClouds.object=struct('Style','checkbox','Value',1,'String','No clouds');
pard.withoutClouds.position=[2,3];
pard.withoutClouds.Width=1;

pard.cloudsThreshold.object=struct('Style','edit','String',0.016);
pard.cloudsThreshold.position=[2,3.6];
pard.cloudsThreshold.Width=0.5;

pard.convertTable.object=struct('Style','text','String','table pos');
pard.convertTable.position=[12,1];
pard.convertTable.Width=2.5;
pard.convertTable.Height=7;

pard.addRule.object=struct('Style','pushbutton','String','+','Callback',{{@addNewRule_callback,obj,'convertTable'}});
pard.addRule.position=[13,1];
pard.addRule.Width=0.2;
pard.addRule.Height=1;

pard.rmRule.object=struct('Style','pushbutton','String','-','Callback',{{@rmRule_callback,obj,'convertTable'}});
pard.rmRule.position=[13,1.2];
pard.rmRule.Width=0.2;
pard.rmRule.Height=1;


% pard.syncParameters={{'roimanager_processors','parsTable',{'Data'}}};

pard.plugininfo.description='Dynamic reconstruction of mammalian CME.';
pard.plugininfo.type='ROI_Analyze';
end
