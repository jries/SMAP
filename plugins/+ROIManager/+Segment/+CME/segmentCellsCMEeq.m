classdef segmentCellsCMEeq<interfaces.DialogProcessor&interfaces.SEProcessor
    %This plug-in is designed for segmentation of yeast cells at their
    %equatorial planes. First, noise will be removed according to their
    %clustered density. Locs will then be converted to pixel images. Next,
    %mask will be generated based on the images and thresholding after
    %gaussian filtering applied. An watershed segmentation algorithm will
    %then be applied. Regions attached to the edge of fields will be
    %further removed. After small objects in the images are removed, cell
    %centers will be set to the mid points of cells along x and y axises.
    methods
        function obj=segmentCellsCMEeq(varargin) 
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={};
        end
        
        function out=run(obj,p)  
          run_segmentCellsCMEeq(obj,p)
          out=[];
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end




function pard=guidef
pard=[];
pard.t1.object=struct('String','Factor of cutoff','Style','text');
pard.t1.position=[1,1];

pard.fOfcutoff.object=struct('String','0.8','Style','edit');
pard.fOfcutoff.position=[1,2];



pard.t2.object=struct('String','Area/pxSize cutoff','Style','text');
pard.t2.position=[2,1];
pard.areaCutoff.object=struct('String','40000','Style','edit');
pard.areaCutoff.position=[2,2];

pard.t3.object=struct('String', 'Density Cutoff', 'Style', 'text');
pard.t3.position=[3,1];
pard.densityCutoff.object = struct('String', '5','Style','edit');
pard.densityCutoff.position=[3,2];

pard.t4.object=struct('String', 'Thresholding cutoff', 'Style', 'text');
pard.t4.position=[4,1];
pard.ampOfIntensityCutOff.object = struct('String', '2','Style','edit');
pard.ampOfIntensityCutOff.position=[4,2];

pard.currentFile.object=struct('String','Current file','Style','checkbox','Value',1);
pard.currentFile.position=[5,1];

pard.preview.object=struct('String','preview','Style','checkbox','Value',1);
pard.preview.position=[6,1];
pard.plugininfo.type='ROI_Segment';

pard.runClusterdensity.object = struct('String', 'Run clusterdensity', 'Style','checkbox', 'Value', 0);
pard.runClusterdensity.position = [7,1];
pard.runClusterdensity.Width = 1;

pard.plugininfo.description='This plug-in is designed for segmentation of yeast cells at their equatorial planes. First, noise will be removed according to their clustered density. Locs will then be converted to pixel images. Next, mask will be generated based on the images and thresholding after gaussian filtering applied. An watershed segmentation algorithm will then be applied. Regions attached to the edge of fields will be further removed. After small objects in the images are removed, cell centers will be set to the mid points of cells along x and y axises.';
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

