classdef ClusterVoronoi<interfaces.DialogProcessor
    %CLUSTERVORONOI Performs Voronoi clustering and returns the density
    %(inverse of voronoi cell). Uses code from SharpViSu:
    %L. Andronov, Y. Lutz, J.-L. Vonesch, and B. P. Klaholz, 
    %?SharpViSu: integrated analysis and segmentation of super-resolution microscopy data,? 
    % Bioinformatics, p. btw123, Mar. 2016.
    methods
        function obj=ClusterVoronoi(varargin)        
            obj@interfaces.DialogProcessor(varargin{:}) ;
            obj.inputParameters={'sr_pixrec','numberOfLayers','sr_pos','sr_size','layers','sr_layerson'};
            obj.history=true;
            obj.showresults=false;
        end
        function out=run(obj,p)
            
            out=vrender(obj,p);
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end
function out=vrender(obj,p)
out=[];
warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId')
[locs,indloc]=obj.locData.getloc({'xnm','ynm','phot','frame','znm','ingrouped','inungrouped'},'layer',find(p.sr_layerson),'position','roi');
A=zeros(length(locs.xnm),9);

A(:,2)=locs.frame;
A(:,4)=locs.xnm-min(locs.xnm);
A(:,5)=locs.ynm-min(locs.ynm);  
A(:,7)=locs.phot;
if ~isempty(locs.znm)
    A(:,6)=locs.znm;
end
Vor = VorArea( A, 1, size(A, 1) );
Area = Vor{1};
ic = Vor{5};
Areacorr = zeros(size(A,1),1);
Areacorr(:) = Area(ic,:);
  
obj.locData.setloc('clusterdensity',Areacorr,indloc);
obj.setPar('locFields',fieldnames(obj.locData));
obj.locData.regroup;
end
           
             
function pard=guidef
pard.t1.object=struct('String','Voronoi cluster analysis from SharpViSu','Style','text');
pard.t1.position=[1,2];
pard.t1.Width=4;

pard.plugininfo.description=sprintf('Performs Voronoi clustering and returns the density (inverse of voronoi cell).If grouped or ungrouped data is used depends on setting in layers. Uses code from SharpViSu: L. Andronov, Y. Lutz, J.-L. Vonesch, and B. P. Klaholz, SharpViSu: integrated analysis and segmentation of super-resolution microscopy data, Bioinformatics, p. btw123, Mar. 2016.');
pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.Name='Cluster Voronoi';
end

