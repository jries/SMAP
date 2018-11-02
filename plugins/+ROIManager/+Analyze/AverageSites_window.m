classdef AverageSites_window<interfaces.DialogProcessor&interfaces.SEProcessor
    methods
        function obj=AverageSites_window(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'se_viewer'};
            obj.showresults=true;
        end
        
        function out=run(obj,p)
            binEdges = p.siteRange(1):p.siteRange(2):p.siteRange(3)+1; 
            roisize=ones(2,1)*obj.P.par.se_siteroi.content;%
            sites=obj.locData.SE.sites;
            
            fdcal=figure(233);
            dcal=plugin('ROIManager','Evaluate','generalStatistics',fdcal,obj.P);
            dcal.attachLocData(obj.SE.locData);
            dcal.makeGui;
            xfirst = roisize(1);
            yfirst = roisize(2);
            
            if p.inOneFile
                sitePerRow = ceil((length(binEdges)-1)/p.rowSites);
                flagPlotCol = 1;
                flagPlotRow = 1;
            end
            newfile = obj.locData.files.filenumberEnd+1;
            oneFilename = ['all_average_sites_F_' num2str(newfile)];
            
            if p.inOneFile&&p.addfile
                obj.locData.addfile(oneFilename);
                initGuiAfterLoad(obj);
                obj.SE.processors.preview.updateFilelist;
            end
                
            for l = 1:length(binEdges)-1
                if sitePerRow<1
                    warning('The row number is larger than the total site number');
                    break
                end
                out=[];
                
                
                if ~p.inOneFile
                    oneFilename = [num2str(binEdges(l)/1000) 'k_' num2str((binEdges(l+1)-1)/1000) 'k'];
                    newfile=obj.locData.files.filenumberEnd+1; % don't update the filenumber if you want to put every ting in
                end
                
                if p.addfile&&~p.inOneFile
                    obj.locData.addfile(oneFilename);
                    initGuiAfterLoad(obj);
                    obj.SE.processors.preview.updateFilelist;
                end
                ticc=tic;
                allIndLoc = [];
                xnmrotAll = [];
                ynmrotAll = [];
                classAll = [];
                for k=binEdges(l):binEdges(l+1)-1 % later to put options here
    %             for k=1:length(sites) % later to put options here
                    if p.sortselection.Value==1 ||  sites(k).annotation.use

                        %PG

                        dcal.site=sites(k);
                        dcal.site.image = obj.locData.SE.plotsite(sites(k));
                        [locsSite,indloc] = dcal.getLocs({'xnmrot','ynmrot'},'size',roisize','grouping','ungrouped');
                        indloc = find(indloc);
                        allIndLoc((end+1):end+length(indloc),:) = indloc;

                        xnmrotAll((end+1):end+length(indloc),:) = locsSite.xnmrot;
                        ynmrotAll((end+1):end+length(indloc),:) = locsSite.ynmrot;
                        classAll((end+1):end+length(indloc),:) = sites(k).ID;
    %                     [locs,indloc]=obj.locData.getloc({'xnm','ynm'},'position',sites(k),'grouping','ungrouped');
        %                 locnew.xnm(indloc)=locnew.xnm(indloc)-sites(k).pos(1);
        %                 locnew.ynm(indloc)=locnew.ynm(indloc)-sites(k).pos(2);
    %                     locnew.xnm(indloc)=locs.xnm-sites(k).pos(1);
    %                     locnew.ynm(indloc)=locs.ynm-sites(k).pos(2);
        %                 figure(88)
        %                 plot(locnew.xnm(indloc),locnew.ynm(indloc),'+')
    %                     locnew.filenumber(indloc)=newfile;
    %                     locnew.class(indloc)=sites(k).ID;
    %                     used=used|indloc;
                    end
                    if toc(ticc)>1
                        ticc=tic;
                        obj.status(['average site: ' num2str(k) ' of ' num2str(length(sites))]); drawnow
                    end
                end
                locnew=obj.locData.loc;
                if p.inOneFile
                    x0=xfirst+(flagPlotCol-1)*roisize(1); % xfirst and yfirst might not be necessary
                    y0=yfirst+(flagPlotRow-1)*roisize(2);
                    flagPlotCol = flagPlotCol+1;
                    if flagPlotCol == sitePerRow+1
                        flagPlotRow = flagPlotRow+1;
                        flagPlotCol = 1;
                    end
                end
                sublocnew = @(x) x(allIndLoc,:);
                locnew = structfun(sublocnew, locnew, 'UniformOutput', false);
                locnew.xnm = xnmrotAll;
                locnew.ynm = ynmrotAll;
                locnew.class = classAll;
                locnew.origin = allIndLoc;
                locnew.filenumber = repelem(newfile, length(allIndLoc))';

                 locc=obj.locData.copy;
                 locc.loc=locnew;
                 locc.regroup;
                 locc.filter;

                if p.addfile
                    locnew.xnm=locnew.xnm+x0;
                    locnew.ynm=locnew.ynm+y0;
                    fn=fieldnames(locnew);
                    for k=1:length(fn)
                        obj.locData.addloc(fn{k},locnew.(fn{k}))
                    end
                    obj.locData.regroup;
                    obj.locData.filter;
                end
           
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
            rr=0:dr:max(max(r1),max(r2))+dr;
            
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

pard.addfile.object=struct('String','add average as new data set','Style','checkbox', 'Value', 1);
pard.addfile.position=[3,3];
pard.addfile.Width=2;

pard.t4.object=struct('String','Range','Style','text');
pard.t4.position=[4,1];
pard.t4.Width=1;

pard.siteRange.object=struct('String','Start binSize end','Style','edit');
pard.siteRange.position=[4,2];
pard.siteRange.Width=1;

pard.t5.object=struct('String','Row of sites','Style','text');
pard.t5.position=[5,1];
pard.t5.Width=1;

pard.rowSites.object=struct('String','1','Style','edit');
pard.rowSites.position=[5,2];
pard.rowSites.Width=1;

pard.inOneFile.object=struct('String','same file','Style','checkbox', 'Value', 1);
pard.inOneFile.position=[5,3];
pard.inOneFile.Width=2;

pard.plugininfo.type='ROI_Analyze';


end