classdef ClusterStatistics_MINFLUX<interfaces.DialogProcessor
%  pair correlation functions calculated according to:
%  Sengupta, Prabuddha, Tijana Jovanovic-Talisman, Dunja Skoko, Malte Renz, 
%  Sarah L Veatch, and Jennifer Lippincott-Schwartz. “Probing Protein Heterogeneity 
%  in the Plasma Membrane Using PALM and Pair Correlation Analysis.” 
%  Nature Methods 8 (September 18, 2011): 969.
    properties
        resultstable
    end
    methods
        function obj=ClusterStatistics_MINFLUX(varargin)    
            obj@interfaces.DialogProcessor(varargin{:}) ;
            obj.inputParameters={'sr_pixrec','numberOfLayers','sr_pos','sr_size','layers','sr_layerson'};
%             obj.history=true;    
            obj.showresults=false;
        end
        
        function out=run(obj,p)
           out=[];          
           if p.filter
                layers=find(obj.getPar('sr_layerson'));
                locs=obj.locData.getloc({'xnm','ynm','znm','clusterindex','tid','groupindex','time','frame','thi'},'layer',layers,'Position','roi');
           else
               % locs=obj.locData.getloc({'xnm','ynm','znm','clusterindex','tid','groupindex','time','frame','thi'},'Position','all','grouping','ungrouped');
               locs=obj.locData.loc;
           end
           

           cfield=p.link.selection;
           if strcmp(cfield,'tid')
               mt=max(locs.tid)+1;
               cfieldv=locs.tid+(single(locs.filenumber)-1)*mt;
           else
               cfieldv=locs.(cfield);
           end
           if isempty(locs.thi)
               dcminflux=false;
           else
               dcminflux=true;
           end

           if isempty(locs.time)
               tfield='frame';
           else
               tfield='time';
           end
           clusterinds=unique(cfieldv(cfieldv>0)); %does not work for multiple files.
           cind=1;
           
           stdsall=zeros(size(obj.locData.loc.xnm),'single');
           stdlall=stdsall;
           angleall=stdsall;
           exts=stdsall;
           extl=stdsall;
           if dcminflux
               previous.time=0;
               previous.cind=0;
               previous.color=0;
               previous.ndc=0;
               previous.indall=[];
               c_dcind=stdsall;
               c_toverlap=stdsall;
               dcind=0;
           end
           for k=1:length(clusterinds)
               inc=cfieldv==clusterinds(k);
               if sum(inc)<p.minloc || sum(inc)-2<=p.skipfirst
                   continue
               end
               indall=cfieldv==clusterinds(k);

               if p.skipfirst>0
                   indo=find(inc,p.skipfirst,'first');
                   inc(indo)=false;
               end

               numlocs(cind,1)=sum(inc);

               xh=double(locs.xnm(inc));yh=double(locs.ynm(inc));

               if dcminflux
                   th=locs.(tfield)(inc);
                   thih=locs.thi(inc);
                   if previous.time(end)>th(1) && thih(1)==1 %overlapping time and other color
                       % if 1 %sufficiently close?
                           if sum(inc)>lencol1 %find longest
                               lencol1=sum(inc); 

                               % time overlap
                               time0=locs.(tfield)(previous.indall);
                               toverlap=min(time0(end),th(end))-max(time0(1),th(1));
                               ttotal=max(time0(end),th(end))-min(time0(1),th(1));
                           % c_dcn1(indall)=sum(inc);
                           % c_dcn1(previous.indall)=previous.ndc;
                               c_toverlap(indall)=toverlap/ttotal;
                               c_toverlap(previous.indall)=toverlap/ttotal;
                               c_dcind(indall)=dcind;
                               c_dcind(previous.indall)=dcind;
                               dcind=dcind+1;
                           end
                       % end
                   end
                   if thih(1)==0 %main color again, go to next
                       previous.time=th;
                       previous.cind=clusterinds(k);
                       previous.color=thih(1);
                       previous.ndc=sum(inc);
                       previous.indall=indall;
                       lencol1=0;
                   end
                   
               end
               xpos(cind,1)=mean(xh);
               ypos(cind,1)=mean(yh);
               
               
               c = cov(xh-mean(xh), yh-mean(yh));
               [a, ev] = eig(c);
               [ev,ind] = sort(diag(ev));
               [xa, ya] = deal(a(1,ind(end)), a(2,ind(end)));
               angle = cart2pol(xa, ya);

               
               stdsall(indall)=sqrt(ev(1));
               stdlall(indall)=sqrt(ev(2));
               angleall(indall)=mod(angle,pi);
               [xr,yr]=rotcoord(xh-mean(xh),yh-mean(yh),angle);
               extl(indall)=max(xr)-min(xr);
               exts(indall)=max(yr)-min(yr);
           end
           if p.addf
               obj.locData.setloc('c_stds',stdsall)
               obj.locData.setloc('c_stdl',stdlall)
               obj.locData.setloc('c_angle',angleall)
               obj.locData.setloc('c_extl',extl)
               obj.locData.setloc('c_exts',exts)
               if dcminflux
                   obj.locData.setloc('c_toverlap',c_toverlap)
                   obj.locData.setloc('c_dcind',c_dcind)
               end

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
pard.link.object=struct('String',{{'clusterindex','groupindex','tid'}},'Style','popupmenu','Value',3);
pard.link.position=[1,1.5];
pard.link.Width=1.5;


pard.filter.object=struct('String','filtered locs (Renderer)','Style','checkbox','Value',0);
pard.filter.position=[1,3];
pard.filter.Width=2;


% pard.advanced.object=struct('String','advanced analysis','Style','checkbox','Value',0);
% pard.advanced.position=[2,1];
% pard.advanced.Width=1.5;

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


pard.plugininfo.name='Cluster Statistics MINFLUX';
pard.plugininfo.description= 'Calculates spatial statistics based on pair correlation and Ripleys K function. pair correlation functions calculated according to: Sengupta, Prabuddha, Tijana Jovanovic-Talisman, Dunja Skoko, Malte Renz, Sarah L Veatch, and Jennifer Lippincott-Schwartz. â€œProbing Protein Heterogeneity  Nature Methods 8 (September 18, 2011): 969.';
pard.plugininfo.type='ProcessorPlugin';

end