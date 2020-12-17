classdef ClusterStatistics<interfaces.DialogProcessor
%  pair correlation functions calculated according to:
%  Sengupta, Prabuddha, Tijana Jovanovic-Talisman, Dunja Skoko, Malte Renz, 
%  Sarah L Veatch, and Jennifer Lippincott-Schwartz. “Probing Protein Heterogeneity 
%  in the Plasma Membrane Using PALM and Pair Correlation Analysis.” 
%  Nature Methods 8 (September 18, 2011): 969.
    properties
        resultstable
    end
    methods
        function obj=ClusterStatistics(varargin)    
            obj@interfaces.DialogProcessor(varargin{:}) ;
            obj.inputParameters={'sr_pixrec','numberOfLayers','sr_pos','sr_size','layers','sr_layerson'};
%             obj.history=true;    
            obj.showresults=true;
        end
        
        function out=run(obj,p)
           out=[];          
           layers=find(obj.getPar('sr_layerson'));
           locs=obj.locData.getloc({'xnm','ynm','znm','clusterindex'},'layer',layers,'Position','roi');
           clusterinds=unique(locs.clusterindex(locs.clusterindex>0));
           cind=1;
           ax=obj.initaxis('boundaries');
           hold(ax,'off')
           plot(ax,locs.xnm,locs.ynm,'b.');
           for k=1:length(clusterinds)
               inc=locs.clusterindex==clusterinds(k);
               if sum(inc)<p.minloc
                   continue
               end
               xh=double(locs.xnm(inc));yh=double(locs.ynm(inc));zh=double(locs.znm(inc));
               numlocs(cind,1)=sum(inc);
               xpos(cind,1)=mean(xh);
               ypos(cind,1)=mean(yh);
               zpos(cind,1)=mean(zh);
               [kc3,volumeCHull(cind,1)]=convhull(xh,yh,zh);
               [kc2,areaCHull(cind,1)]=convhull(xh,yh);
               [kb3,volumeB(cind,1)]=boundary(xh,yh,zh);
               [kb2,areaB(cind,1)]=boundary(xh,yh);
               cID(cind,1)=clusterinds(k);
               cind=cind+1;
               hold(ax,'on')
               plot(ax,xh,yh,'+')
               
               plot(ax,xh(kc2),yh(kc2),'r')
               plot(ax,xh(kb2),yh(kb2),'m')
           end
           obj.resultstable=table(cID,numlocs,areaCHull,volumeCHull,areaB,volumeB,xpos,ypos,zpos);
           
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function exportbutton(obj,a,b)
            tableo=obj.resultstable;
            fn=obj.getPar('lastSMLFile');
            if ~isempty(fn)
                fno=strrep(fn,'_sml.mat','_cluster.csv');
            else
                fno='*.csv';
            end
           [f,p]=uiputfile(fno);
           if f
               writetable(tableo,[p f])
           end
            
        end

    end
end




function pard=guidef(obj)
pard.minloct.object=struct('String','min number of locs per cluster','Style','text');
pard.minloct.position=[1,1];
pard.minloct.Width=2;
pard.minloc.object=struct('String','10','Style','edit');
pard.minloc.position=[1,3];
pard.minloc.Width=0.5;

pard.export.object=struct('String','export','Style','pushbutton','Callback',@obj.exportbutton);
pard.export.position=[2,1];
pard.export.Width=1;


pard.plugininfo.name='Cluster Statistics';
pard.plugininfo.description= 'Calculates spatial statistics based on pair correlation and Ripleys K function. pair correlation functions calculated according to: Sengupta, Prabuddha, Tijana Jovanovic-Talisman, Dunja Skoko, Malte Renz, Sarah L Veatch, and Jennifer Lippincott-Schwartz. â€œProbing Protein Heterogeneity  Nature Methods 8 (September 18, 2011): 969.';
pard.plugininfo.type='ProcessorPlugin';

end