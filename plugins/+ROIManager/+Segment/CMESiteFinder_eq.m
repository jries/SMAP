classdef CMESiteFinder_eq<interfaces.DialogProcessor&interfaces.SEProcessor
    methods
        function obj=CMESiteFinder_eq(varargin) 
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={};
        end
        
        function out=run(obj,p)  
%           run_CMESiteFinder_eq(obj,p)
          run_CMESiteFinder_eq_new(obj,p)
          out=[];
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end




function pard=guidef
pard=[];
pard.t_width_widerZone.object=struct('String','Wide zone width (nm)','Style','text');
pard.t_width_widerZone.position=[1,1];
pard.width_widerZone.object=struct('String','250','Style','edit');
pard.width_widerZone.position=[1,2];

pard.t_width_narrowerZone.object=struct('String','Narrow zone width (nm)','Style','text');
pard.t_width_narrowerZone.position=[2,1];
pard.width_narrowerZone.object=struct('String','70','Style','edit');
pard.width_narrowerZone.position=[2,2];

pard.t_minNumOfLocs_siteRemoval.object=struct('String','Min # of locs per site','Style','text');
pard.t_minNumOfLocs_siteRemoval.position=[3,1];
pard.minNumOfLocs_siteRemoval.object=struct('String','20','Style','edit');
pard.minNumOfLocs_siteRemoval.position=[3,2];

pard.t_epsilon_segmentation.object=struct('String','DBSCAN Epsilon','Style','text');
pard.t_epsilon_segmentation.position=[1,3];
pard.epsilon_segmentation.object=struct('String','50','Style','edit');
pard.epsilon_segmentation.position=[1,4];

pard.t_MinPts_segmentation.object=struct('String','DBSCAN min density','Style','text');
pard.t_MinPts_segmentation.position=[2,3];
pard.MinPts_segmentation.object=struct('String','9','Style','edit');
pard.MinPts_segmentation.position=[2,4];

pard.t_dist_fromBound.object=struct('String','Offset from boundary (nm)','Style','text');
pard.t_dist_fromBound.position=[3,3];
pard.dist_fromBound.object=struct('String','70','Style','edit');
pard.dist_fromBound.position=[3,4];

pard.t_length_tangentLine.object=struct('String','Tengent line length (nm)','Style','text');
pard.t_length_tangentLine.position=[4,3];
pard.length_tangentLine.object=struct('String','200','Style','edit');
pard.length_tangentLine.position=[4,4];

pard.preview.object=struct('String','preview selected cell','Style','checkbox','Value',1);
pard.preview.position=[6,1];

pard.plugininfo.type='ROI_Analyze';


end

