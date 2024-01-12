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
           locs=obj.locData.getloc({'xnm','ynm','znm','clusterindex','tid','groupindex'},'layer',layers,'Position','roi');
           cfield=p.link.selection;
           clusterinds=unique(locs.(cfield)(locs.(cfield)>0));
           cind=1;
           if p.advanced
               ax=obj.initaxis('boundaries');
               hold(ax,'off')
               plot(ax,locs.xnm,locs.ynm,'b.');
           end
           stdsall=zeros(size(obj.locData.loc.xnm),'single');
           stdlall=stdsall;
           angleall=stdsall;
           exts=stdsall;
           extl=stdsall;
           for k=1:length(clusterinds)
               inc=locs.(cfield)==clusterinds(k);
               if sum(inc)<p.minloc || sum(inc)-2<=p.skipfirst
                   continue
               end
               

               if p.skipfirst>0
                   indo=find(inc,p.skipfirst,'first');
                   inc(indo)=false;
               end

               numlocs(cind,1)=sum(inc);

               xh=double(locs.xnm(inc));yh=double(locs.ynm(inc));
               
               xpos(cind,1)=mean(xh);
               ypos(cind,1)=mean(yh);
               
               if p.advanced
                   if ~isempty(locs.znm)
                       zh=double(locs.znm(inc));
                       zpos(cind,1)=mean(zh);
                       [kc3,volumeCHull(cind,1)]=convhull(xh,yh,zh);
                       [kb3,volumeB(cind,1)]=boundary(xh,yh,zh);
                   else
                       volumeCHull(cind,1)=0;
                       volumeB(cind,1)=0;
                       zpos(cind,1)=0;
                   end
                   [kc2,areaCHull(cind,1)]=convhull(xh,yh);
                   [kb2,areaB(cind,1)]=boundary(xh,yh);
                   cID(cind,1)=clusterinds(k);
                   cind=cind+1;
                   hold(ax,'on')
                   plot(ax,xh,yh,'+')
                   
                   plot(ax,xh(kc2),yh(kc2),'r')
                   plot(ax,xh(kb2),yh(kb2),'m')
               end
                
               
               c = cov(xh-mean(xh), yh-mean(yh));
               [a, ev] = eig(c);
               [ev,ind] = sort(diag(ev));
               [xa, ya] = deal(a(1,ind(end)), a(2,ind(end)));
               angle = cart2pol(xa, ya);

               indall=obj.locData.loc.(cfield)==clusterinds(k);
               stdsall(indall)=sqrt(ev(1));
               stdlall(indall)=sqrt(ev(2));
               angleall(indall)=mod(angle,pi);
               [xr,yr]=rotcoord(xh-mean(xh),yh-mean(yh),angle);
               extl(indall)=max(xr)-min(xr);
               exts(indall)=max(yr)-min(yr);
           end
           if p.advanced
            obj.resultstable=table(cID,numlocs,areaCHull,volumeCHull,areaB,volumeB,xpos,ypos,zpos);
           else
           end
           if p.addf
               obj.locData.setloc('c_stds',stdsall)
               obj.locData.setloc('c_stdl',stdlall)
               obj.locData.setloc('c_angle',angleall)
               obj.locData.setloc('c_extl',extl)
               obj.locData.setloc('c_exts',exts)
               obj.locData.regroup;
           end

           
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
pard.minloct.position=[3,1];
pard.minloct.Width=2;
pard.minloc.object=struct('String','10','Style','edit');
pard.minloc.position=[3,2.5];
pard.minloc.Width=0.5;


pard.linkt.object=struct('String','Selection','Style','text');
pard.linkt.position=[1,1];
pard.link.object=struct('String',{{'clusterindex','groupindex','tid'}},'Style','popupmenu','Value',2);
pard.link.position=[1,2];
pard.link.Width=1.5;

pard.advanced.object=struct('String','advanced analysis','Style','checkbox','Value',1);
pard.advanced.position=[2,1];
pard.advanced.Width=1.5;

pard.addf.object=struct('String','add results to locs','Style','checkbox','Value',0);
pard.addf.position=[2,3];
pard.addf.Width=1.5;

pard.skipt.object=struct('String','Skip first','Style','text');
pard.skipt.position=[3,3];
pard.skipfirst.object=struct('String',0,'Style','edit');
pard.skipfirst.position=[3,4];
pard.skipfirst.Width=0.5;



pard.export.object=struct('String','export','Style','pushbutton','Callback',@obj.exportbutton);
pard.export.position=[5,1];
pard.export.Width=1;


pard.plugininfo.name='Cluster Statistics';
pard.plugininfo.description= 'Calculates spatial statistics based on pair correlation and Ripleys K function. pair correlation functions calculated according to: Sengupta, Prabuddha, Tijana Jovanovic-Talisman, Dunja Skoko, Malte Renz, Sarah L Veatch, and Jennifer Lippincott-Schwartz. â€œProbing Protein Heterogeneity  Nature Methods 8 (September 18, 2011): 969.';
pard.plugininfo.type='ProcessorPlugin';

end