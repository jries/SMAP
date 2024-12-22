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
            obj.showresults=true;
        end
        
        function out=run(obj,p)
           out=[];          
           if p.filter
                layers=find(obj.getPar('sr_layerson'));
                locs=obj.locData.getloc({'xnm','ynm','znm','clusterindex','tid','groupindex','time','frame','thi','filenumber','phot'},'layer',layers,'Position','roi');
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

           % statistics per cluster cstat.xxx. 
           maxid=max(clusterinds);
           cstat.numlocs=zeros(maxid,1);
           cstat.stdx=zeros(maxid,1);
           cstat.stdy=zeros(maxid,1);
           cstat.xpos=zeros(maxid,1);
           cstat.ypos=zeros(maxid,1);
           cstat.binx=zeros(maxid,15);
           cstat.biny=zeros(maxid,15);
           cstat.taillength=zeros(maxid,1);
           cstat.extl=zeros(maxid,1);
           cstat.exts=zeros(maxid,1);
           cstat.stdlong=zeros(maxid,1);
           cstat.stdshort=zeros(maxid,1);
           cstat.timediffmedian=zeros(maxid,1);
           cstat.binrelx=zeros(maxid,15);
           cstat.binrely=zeros(maxid,15); 
           cstat.binphot=zeros(maxid,15);
           cstat.sigminflux=zeros(maxid,15);

           for k=1:length(clusterinds)
               cind=clusterinds(k);
               indall=cfieldv==cind; %all track
               inc=indall; % 
               if sum(indall)<p.minloc || sum(indall)-2<=p.skipfirst
                   continue
               end

               if p.skipfirst>0
                   indo=find(inc,p.skipfirst,'first');
                   inc(indo)=false;
               end

               cstat.numlocs(cind,1)=sum(inc);

               xh=double(locs.xnm(inc));yh=double(locs.ynm(inc));
               photh=double(locs.phot(inc));

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
               th=double(locs.time(inc));
               dt=diff(th);
               cstat.timediffmedian(cind)=median(dt);


               cstat.xpos(cind,1)=mean(xh);
               cstat.ypos(cind,1)=mean(yh);
               
               
               c = cov(xh-mean(xh), yh-mean(yh));
               [a, ev] = eig(c);
               [ev,ind] = sort(diag(ev));
               [xa, ya] = deal(a(1,ind(end)), a(2,ind(end)));
               angle = cart2pol(xa, ya);

               stds=sqrt(ev(1));
               stdl=sqrt(ev(2));

               stdsall(indall)=stds;
               stdlall(indall)=stdl;
               angleall(indall)=mod(angle,pi);
               [xr,yr]=rotcoord(xh-mean(xh),yh-mean(yh),angle);
               cstat.extl(cind)=max(xr)-min(xr);
               cstat.exts(cind)=max(yr)-min(yr);

               cstat.stdx(cind)=std(xh);
               cstat.stdy(cind)=std(yh);
               cstat.stdlong(cind)=stdl;
               cstat.stdshort(cind)=stds;

               %binning
               L=75;%?
               [stdx,photb,sigminflux]=bintraceMINFLUX(xh, photh, L);
               sxrel=stdx./sigminflux;
               cstat.binx(cind,1:min(15,length(stdx)))=stdx(1:min(15,length(stdx)));
               cstat.binrelx(cind,1:min(15,length(stdx)))=sxrel(1:min(15,length(stdx)));
               [stdy,~]=bintraceMINFLUX(yh, photh, L);
               syrel=stdy./sigminflux;
               cstat.biny(cind,1:min(15,length(stdy)))=stdy(1:min(15,length(stdy)));
               cstat.binrely(cind,1:min(15,length(stdy)))=syrel(1:min(15,length(stdy)));
               cstat.binphot(cind,1:min(15,length(stdy)))=photb(1:min(15,length(stdy)));
               cstat.sigminflux(cind,1:min(15,length(sigminflux)))=sigminflux(1:min(15,length(sigminflux)));

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

           %plotting
           axstdx=obj.initaxis('stdx');
           stdxp=cstat.stdx(cstat.numlocs>0);
           stdyp=cstat.stdy(cstat.numlocs>0);
           numlocs=cstat.numlocs(cstat.numlocs>0);
           
           nhist=0:1:max(quantile(stdxp,0.95),quantile(stdyp,0.95));
           
           ff='%2.1f';
           histogram(axstdx,stdxp,nhist);
           title(axstdx,['median ' num2str(median(stdxp),ff)])
           axstdx.XLim(1)=0;
           xlabel(axstdx,'std(x) nm')

           axstdy=obj.initaxis('stdy');
           histogram(axstdy,stdyp,nhist);
           title(axstdy,['median ' num2str(median(stdyp),ff)])
           axstdy.XLim(1)=0;
           xlabel(axstdy,'std(y) nm')

           axstdy=obj.initaxis('locs');
           histogram(axstdy,numlocs);
           title(axstdy,['median ' num2str(median(numlocs),ff)])
           axstdy.XLim(1)=0;
           xlabel(axstdy,'localizations / track')

           axdt=obj.initaxis('dt');
           histogram(axdt,cstat.timediffmedian(cstat.numlocs>0));
           title(axdt,['median ' num2str(median(cstat.timediffmedian(cstat.numlocs>0)),ff)])
           axdt.XLim(1)=0;
           xlabel(axdt,'median time diff (ms)')


           quantvals=[0.25, 0.5, 0.75];
           for k=size(cstat.binrelx,2):-1:1
               v=cstat.binrelx(:,k);
               if sum(v>0)<5
                   continue
               end
               temp=quantile(v(v>0),quantvals,1);
               binxmed(k)=temp(2); binxrel05(k)=temp(1); binxrel95(k)=temp(3);
               v=cstat.binrely(:,k);
               temp=quantile(v(v>0),quantvals,1);
               binymed(k)=temp(2); binyrel05(k)=temp(1); binyrel95(k)=temp(3);
               v=cstat.binphot(:,k);
               temp=quantile(v(v>0),quantvals,1);
               binphotmed(k)=temp(2);

               v=cstat.binx(:,k);
               temp=quantile(v(v>0),quantvals,1);
               binx(k)=temp(2); binx05(k)=temp(1); binx95(k)=temp(3);
               v=cstat.biny(:,k);
               temp=quantile(v(v>0),quantvals,1);
               biny(k)=temp(2); biny05(k)=temp(1); biny95(k)=temp(3);
               v=cstat.sigminflux(:,k);
               temp=quantile(v(v>0),quantvals,1);
               sigminfluxmed(k)=temp(2); 

           end
           axbx=obj.initaxis('binrel');
           hold(axbx,'off')
           semilogx(axbx,binphotmed,binxmed,'b+-',binphotmed,binxrel05,'b--',binphotmed,binxrel95,'b--')
           hold(axbx,'on')
           semilogx(axbx,binphotmed,binymed,'r+-',binphotmed,binyrel05,'r--',binphotmed,binyrel95,'r--')
           xlabel(axbx,'median photons')
           ylabel(axbx,'std/locprec')
           legend(axbx,'x: median',[num2str(quantvals(1)*100) '%'],[num2str(quantvals(3)*100) '%'],'y median')

           axbin=obj.initaxis('bin');
           hold(axbin,'off')
           semilogx(axbin,binphotmed,binx,'b+-',binphotmed,binx05,'b--',binphotmed,binx95,'b--')
           hold(axbin,'on')
           semilogx(axbin,binphotmed,sigminfluxmed,'k')
           semilogx(axbin,binphotmed,biny,'r+-',binphotmed,biny05,'r--',binphotmed,biny95,'r--')
           


           xlabel(axbin,'median photons')
           ylabel(axbin,'std (nm)')
           legend(axbin,'x: median',[num2str(quantvals(1)*100) '%'],[num2str(quantvals(3)*100) '%'],'locprec MF', 'y median')
           
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