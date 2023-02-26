classdef segmentCellsCMEeq<interfaces.DialogProcessor&interfaces.SEProcessor
% This plug-in is designed for segmentation of yeast cells at their
% equatorial planes. 
    methods
        function obj=segmentCellsCMEeq(varargin) 
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={};
        end
        
        function out=run(obj,p)
            run_segmentCellsCMEeq_new(obj,p)
          out=[];
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end


function pard=guidef
pard=[];

pard.t_win_partition.object=struct('String','Partition size (nm)','Style','text');
pard.t_win_partition.position=[1,1];
pard.win_partition.object=struct('String','1000','Style','edit');
pard.win_partition.position=[1,2];

pard.t_margin_partition.object=struct('String', 'Partition overlap (nm)', 'Style', 'text');
pard.t_margin_partition.position=[2,1];
pard.margin_partition.object = struct('String', '50','Style','edit');
pard.margin_partition.position=[2,2];

pard.t_epsilon_segmentation.object=struct('String', 'DBSCAN Epsilon', 'Style', 'text');
pard.t_epsilon_segmentation.position=[3,1];
pard.epsilon_segmentation.object = struct('String', '120','Style','edit');
pard.epsilon_segmentation.position=[3,2];

pard.t_MinPts_segmentation.object=struct('String', 'DBSCAN min density', 'Style', 'text');
pard.t_MinPts_segmentation.position=[4,1];
pard.MinPts_segmentation.object = struct('String', '5','Style','edit');
pard.MinPts_segmentation.position=[4,2];

pard.t_edgeMargin_cellRemoval.object=struct('String', 'Fov padding (nm)', 'Style', 'text');
pard.t_edgeMargin_cellRemoval.position=[1,3];
pard.edgeMargin_cellRemoval.object = struct('String', '10','Style','edit');
pard.edgeMargin_cellRemoval.position=[1,4];

pard.t_minNumOfLocs_cellRemoval.object=struct('String', 'Min # of locs per cell', 'Style', 'text');
pard.t_minNumOfLocs_cellRemoval.position=[2,3];
pard.minNumOfLocs_cellRemoval.object = struct('String', '1000','Style','edit');
pard.minNumOfLocs_cellRemoval.position=[2,4];

pard.t_zoneWidth_boundary.object=struct('String', 'Boundary zone width (nm)', 'Style', 'text');
pard.t_zoneWidth_boundary.position=[3,3];
pard.zoneWidth_boundary.object = struct('String', '30','Style','edit');
pard.zoneWidth_boundary.position=[3,4];

pard.t_shrinkFactor_boundary.object=struct('String', 'Boundary shrink factor', 'Style', 'text');
pard.t_shrinkFactor_boundary.position=[4,3];
pard.shrinkFactor_boundary.object = struct('String', '0.98','Style','edit');
pard.shrinkFactor_boundary.position=[4,4];

pard.t_smoothFactor_boundary.object=struct('String', 'Boundary smooth factor', 'Style', 'text');
pard.t_smoothFactor_boundary.position=[5,3];
pard.smoothFactor_boundary.object = struct('String', '0.99994','Style','edit');
pard.smoothFactor_boundary.position=[5,4];

pard.preview.object=struct('String','Preview selected file','Style','checkbox','Value',1);
pard.preview.position=[7,1];
pard.preview.Tooltip = 'Preview the cell segmentation applied to the currecnt file selected in the ROI manager.';

pard.plugininfo.type='ROI_Segment';
pard.plugininfo.description='This plug-in is designed for segmentation of yeast cells at their equatorial planes.';
% 
% pard.t2.object=struct('String','sigmaNMS','Style','text');
% pard.t2.position=[2,1];
% pard.sigmaNMS.object=struct('String','5','Style','edit');
% pard.sigmaNMS.position=[2,2];
% pard.t3.object=struct('String','diameterNPC','Style','text');
% pard.t3.position=[3,1];
% pard.diameterNPC.object=struct('String','110','Style','edit');
% pard.diameterNPC.position=[3,2];
% pard.t4.object=struct('String','rim','Style','text');
% pard.t4.position=[4,1];
% pard.rim.object=struct('String','20','Style','edit');
% pard.rim.position=[4,2];
% 
% pard.saveon.object=struct('String','saveon','Style','checkbox');
% pard.saveon.position=[1,3];
% 
% pard.getmask.object=struct('String','getmask','Style','checkbox');
% pard.getmask.position=[2,3];
end

