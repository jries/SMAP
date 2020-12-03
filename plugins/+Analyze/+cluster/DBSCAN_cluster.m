classdef DBSCAN_cluster<interfaces.DialogProcessor
    %   Copyright (c)2017 Ries Lab, European Molecular Biology Laboratory,
    %   Heidelberg. This file is part of Single Molecule Analysis Platform (SMAP).
    %   DBSCAN_cluster is an optimized version of the DBSCAN algorithm from 
    %   F. A. Caetano, B. S. Dirk, J. H. K. Tam, P. C. Cavanagh, M. Goiko, S. S. G. Ferguson, S. H. Pasternak, 
    %   J. D. Dikeakos, J. R. de Bruyn, and B. Heit, ?MIiSR: Molecular Interactions in Super-Resolution Imaging 
    %   Enables the Analysis of Protein Interactions, Dynamics and Formation of Multi-protein Structures.,? 
    %   PLoS Comput. Biol., vol. 11, no. 12, p. e1004634, Dec. 2015.
    methods
        function obj=DBSCAN_cluster(varargin)        
            obj@interfaces.DialogProcessor(varargin{:}) ;
            obj.inputParameters={'sr_pixrec','numberOfLayers','sr_pos','sr_size','layers','sr_layerson'};
            obj.history=false;
            obj.showresults=true;
        end
        function out=run(obj,p)
            [locs,indused]=obj.locData.getloc({'xnm','ynm','znm'},'layer',find(p.sr_layerson),'position','roi');
            if p.seteps
                eps=p.eps_dbscan;
            else
                eps=[];
            end
            setstatus=@obj.status;
            ax=obj.initaxis('DBSCAN');
            if isempty(locs.znm)
                locs=rmfield(locs,'znm');
            end
            [DBSCANmat,DBSCANtab]=DBSCANsmap(locs,p.k_dbscan,eps,ax,setstatus);
            clusterind=DBSCANtab(:,4);
            obj.locData.setloc('clusterindex',clusterind,indused);
            obj.locData.setloc('clustertype',DBSCANtab(:,5),indused);
            obj.setPar('locFields',fieldnames(obj.locData.loc));
            out=DBSCANmat;
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end
 
             
function pard=guidef
pard.t1.object=struct('String','DBSCAN cluster analysis from MiISR','Style','text');
pard.t1.position=[1,1];
pard.t1.Width=4;

pard.t2.object=struct('String','minimum objects in neighoburhood (k)','Style','text');
pard.t2.position=[2,1];
pard.t2.Width=3;

pard.k_dbscan.object=struct('String','10','Style','edit');
pard.k_dbscan.position=[2,4];

pard.seteps.object=struct('String','set eps (neighbourhood radius) to: ','Style','checkbox','Value',0);
pard.seteps.position=[3,1];
pard.seteps.Width=3;

pard.eps_dbscan.object=struct('String','5','Style','edit');
pard.eps_dbscan.position=[3,4];

pard.plugininfo.description=sprintf(['DBSCAN_cluster is an optimized version of the DBSCAN algorithm from',...
    'F. A. Caetano, B. S. Dirk, J. H. K. Tam, P. C. Cavanagh, M. Goiko, S. S. G. Ferguson, S. H. Pasternak,',... 
    'J. D. Dikeakos, J. R. de Bruyn, and B. Heit, ?MIiSR: Molecular Interactions in Super-Resolution Imaging ',...
    'Enables the Analysis of Protein Interactions, Dynamics and Formation of Multi-protein Structures.,? ',...
    'PLoS Comput. Biol., vol. 11, no. 12, p. e1004634, Dec. 2015.'    ]);
pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.name='Cluster DBSCAN';

end

