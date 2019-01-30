classdef CMESiteFinder_eq<interfaces.DialogProcessor&interfaces.SEProcessor
    methods
        function obj=CMESiteFinder_eq(varargin) 
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={};
        end
        
        function out=run(obj,p)  
          run_CMESiteFinder_eq(obj,p)
          out=[];
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end




function pard=guidef
pard=[];
pard.t1.object=struct('String','Density cut off','Style','text');
pard.t1.position=[1,1];
pard.th.object=struct('String','0.55','Style','edit');
pard.th.position=[1,2];

pard.t2.object=struct('String','Bandwidth','Style','text');
pard.t2.position=[2,1];
pard.bw.object=struct('String','0.01','Style','edit');
pard.bw.position=[2,2];

pard.t3.object=struct('String','Circularity','Style', 'text');
pard.t3.position=[3,1];

pard.filter_circularity.object = struct('String', '', 'Style', 'checkbox', 'Value', 1);
pard.filter_circularity.position = [3,1.8];

pard.cutoff_circularity.object = struct('String', '0.9', 'Style', 'edit');
pard.cutoff_circularity.position = [3,2];

pard.preview.object=struct('String','preview','Style','checkbox','Value',1);
pard.preview.position=[8,1];
pard.plugininfo.type='ROI_Analyze';
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

