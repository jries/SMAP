classdef AverageSites<interfaces.DialogProcessor&interfaces.SEProcessor
%     Calculates an average of selected ROIs by overlaying them and
%     plotting all localizations into a new file.
    methods
        function obj=AverageSites(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'se_viewer'};
            obj.showresults=true;
        end
        
        function out=run(obj,p)  
            out=[];
            
            if p.addfile
                obj.locData.addfile(p.name);
                initGuiAfterLoad(obj);
                obj.SE.processors.preview.updateFilelist;
            end
            newfile=obj.locData.files.filenumberEnd;
            locnew=obj.locData.loc;
            sites=obj.locData.SE.sites;
            used=false(size(locnew.xnm));
            x0=nanmedian(locnew.xnm);
            y0=nanmedian(locnew.ynm);
            if ~isfield(locnew,'class')
                locnew.class=0*locnew.xnm;
            end
            ticc=tic;
            for k=1:length(sites)
                if p.sortselection.Value==1 ||  sites(k).annotation.use
                    [locs,indloc]=obj.locData.getloc({'xnm','ynm'},'position',sites(k),'grouping','ungrouped');
    %                 locnew.xnm(indloc)=locnew.xnm(indloc)-sites(k).pos(1);
    %                 locnew.ynm(indloc)=locnew.ynm(indloc)-sites(k).pos(2);
                    locnew.xnm(indloc)=locs.xnm-sites(k).pos(1);
                    locnew.ynm(indloc)=locs.ynm-sites(k).pos(2);
    %                 figure(88)
    %                 plot(locnew.xnm(indloc),locnew.ynm(indloc),'+')
                    locnew.filenumber(indloc)=newfile;
                    locnew.class(indloc)=sites(k).ID;
                    used=used|indloc;
                end
                if toc(ticc)>1
                    ticc=tic;
                    obj.status(['average site: ' num2str(k) ' of ' num2str(length(sites))]); drawnow
                end
            end
            
            used=used& ~( isnan(locnew.xnm) | isnan(locnew.ynm) | isnan(locnew.class));
            
            fn=fieldnames(locnew);
              for k=1:length(fn)
                   locnew.(fn{k})=locnew.(fn{k})(used);
              end
             
             locc=obj.locData.copy;
             locc.loc=locnew;
             locc.regroup;
             locc.filter;

            if p.addfile
                locnew.xnm=locnew.xnm+x0;
                locnew.ynm=locnew.ynm+y0;
                obj.locData.addLocData(locnew);
%                 for k=1:length(fn)
%                     obj.locData.addloc(fn{k},locnew.(fn{k}))
%                 end
                obj.locData.regroup;
                obj.locData.filter;
            end
           
            
            
            %try: add empty file, there put averaged sites
            %for every site: loc.xnm-site.pos(1)+xpossite
            
            loc1=locc.getloc({'xnm','ynm'},'layer',1,'position','all','filenumber',newfile);
            loc2=locc.getloc({'xnm','ynm'},'layer',2,'position','all','filenumber',newfile);
            %do some statistics:
           
            
            
            ax=obj.initaxis('scatter');
            plot(loc1.xnm,loc1.ynm,'.')
            hold on
            plot(loc2.xnm,loc2.ynm,'.')
            hold off
            
             [phi1,r1]=cart2pol(loc1.xnm,loc1.ynm);
             [phi2,r2]=cart2pol(loc2.xnm,loc2.ynm);
            ax=obj.initaxis('radial distribution');
            dr=5;
%             rr=dr/2:dr:max(max(r1),max(r2));
            %proper concentration: dA=pi*(ro^2-ri^2)
            rr=0:dr:max(max([r1;r2]))+dr;
            
            dA=pi*(rr(2:end)+rr(1:end-1)).*(rr(2:end)-rr(1:end-1));
            h1=histcounts(r1,[rr]) ;
             h2=histcounts(r2,[rr]) ;
              rrp=rr(1:end-1)+dr/2;
            plot(rrp,h1,rrp,h2)
            xlabel('r')
            ylabel('counts')
            
            ax=obj.initaxis('radial concentration');
%             plot(rr,h1./rr/2/pi,rr,h2./rr/2/pi)
           
             plot(rrp,h1./dA,rrp,h2./dA)
            xlabel('r')
            ylabel('concentration')
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end




function pard=guidef

pard.t1.object=struct('String','average sites','Style','text');
pard.t1.position=[1,1];
pard.t1.Width=1;

pard.post.object=struct('String','all or annotated:use','Style','text');
pard.post.position=[2,1];
pard.post.Width=1;

% pard.pos.object=struct('String','0,0','Style','edit');
% pard.pos.position=[2,2];
% pard.pos.Width=1;


pard.sortselection.object=struct('String',{{'all','use'}},'Style','popupmenu');
pard.sortselection.position=[2,2];
pard.sortselection.Width=1;

pard.namet.object=struct('String','name','Style','text');
pard.namet.position=[3,1];
pard.namet.Width=1;

pard.name.object=struct('String','average','Style','edit');
pard.name.position=[3,2];
pard.name.Width=1;

pard.addfile.object=struct('String','add average as new data set','Style','checkbox');
pard.addfile.position=[3,3];
pard.addfile.Width=2;

pard.plugininfo.type='ROI_Analyze';
pard.plugininfo.description='Calculates an average of selected ROIs by overlaying them and plotting all localizations into a new file.';

end