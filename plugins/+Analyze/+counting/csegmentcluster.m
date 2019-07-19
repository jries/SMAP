classdef csegmentcluster<interfaces.DialogProcessor
%     Watershed-based segmentation of clusters
    methods
        function obj=csegmentcluster(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});  
        end
        
        function out=run(obj,p)           
            locs=obj.locData.getloc({'frame','xnm','ynm','phot','bg','PSFxnm','locprecnm'},'layer',1,'position','roi');

%             clusters=cluster_counting(locs,p);
            clusters=cluster_counting_2(obj.locData,p);
            obj.setResults('counting_clusters',clusters);
            out=clusters;
        end
        
        function refit_callback(obj)
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end




function pard=guidef(obj)
% 
pard.text2.object=struct('String','length scale nm','Style','text');
pard.text2.position=[1,1];

pard.lengthscale.object=struct('String','15','Style','edit');
pard.lengthscale.position=[1,2];
pard.lengthscale.isnumeric=1;

pard.text3.object=struct('String','maxPSF in focus','Style','text');
pard.text3.position=[2,1];
pard.segment_maxPSF.object=struct('String','130','Style','edit');
pard.segment_maxPSF.position=[2,2];
pard.segment_maxPSF.isnumeric=1;

pard.plugininfo.name='segment cluster';
pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.description='Watershed-based segmentation of clusters';
end