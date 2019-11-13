classdef AverageSites_NPC<interfaces.DialogProcessor&interfaces.SEProcessor
    methods
        function obj=AverageSites_NPC(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'se_viewer'};
            obj.showresults=true;
        end
        
        function out=run(obj,p)
            binEdges = p.siteRange(1):p.siteRange(2):p.siteRange(3)+1; % = how many rois/bin
            roisize=ones(2,1)*obj.P.par.se_siteroi.content+200;
            sites=obj.locData.SE.sites;
            
            % hack an evaluate plug-in in order to use the obj.getLocs(...)
            fdcal=figure(233);
            dcal=plugin('ROIManager','Evaluate','generalStatistics',fdcal,obj.P);
            dcal.attachLocData(obj.SE.locData);
            dcal.makeGui;
            xfirst = roisize(1);
            yfirst = roisize(2);
            
            if p.inOneFile
            % called when the use wants to put all average ROIs into the same file
                sitePerRow = ceil((length(binEdges)-1)/p.rowSites);
                flagPlotCol = 1;
                flagPlotRow = 1;
            end
            newfile = obj.locData.files.filenumberEnd+1;
            oneFilename = ['all_average_sites_F_' num2str(newfile)];
            
            if p.inOneFile&&p.addfile
             % called when the use wants to save and to put all average ROIs into the same file
                obj.locData.addfile(oneFilename);
                initGuiAfterLoad(obj);
                obj.SE.processors.preview.updateFilelist;
            end
                
            for l = 1:length(binEdges)-1
                % per bin execution
                if sitePerRow<1
                    warning('The row number is larger than the total site number');
                    break
                end
                out=[];
                
                
                if ~p.inOneFile
                    oneFilename = [num2str(binEdges(l)/1000) 'k_' num2str((binEdges(l+1)-1)/1000) 'k'];
                    newfile=obj.locData.files.filenumberEnd+1; % update the filenumber/new file
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
                orderAll = [];
                znmAll = [];
                channelAll = [];
                for k=binEdges(l):binEdges(l+1)-1 % later to put options here
                  % per site execution
    %             for k=1:length(sites) % later to put options here
                    if p.sortselection.Value==1 ||  sites(k).annotation.use

                        % Borrow the evaluate plug-in to use getLocs(obj,...)
                        dcal.site=sites(k);
                        dcal.site.image = obj.locData.SE.plotsite(sites(k));
                        [locsSite,indloc] = dcal.getLocs({'xnmrot','ynmrot','znm','channel'},'size',roisize','grouping','ungrouped'); % per ROI info.
                        fitting = obj.getPar(['fitter_' 'SMLMModelFitGUI2']);
                        fitting.allParsArg = sites(k).evaluation.SMLMModelFitGUI2.allParsArg;
                        locsSite.xnm = locsSite.xnmrot;
                        locsSite.ynm = locsSite.ynmrot;
                        locsSiteNew = fitting.locsHandler(locsSite, fitting.exportPars(1,'lPar'));
                        indloc = find(indloc); 
                        allIndLoc((end+1):end+length(indloc),:) = indloc; % all ROIs Ind.

                        xnmrotAll((end+1):end+length(indloc),:) = locsSiteNew.xnm;
                        ynmrotAll((end+1):end+length(indloc),:) = locsSiteNew.ynm;
                        znmAll((end+1):end+length(indloc),:) = locsSiteNew.znm;
                        channelAll((end+1):end+length(indloc),:) = locsSiteNew.channel;
                        classAll((end+1):end+length(indloc),:) = sites(k).ID;
                        orderAll((end+1):end+length(indloc),:) = k;
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
                sublocnew = @(x) x(allIndLoc,:); % subset the sites
                locnew = structfun(sublocnew, locnew, 'UniformOutput', false);
                locnew.xnm = xnmrotAll; % updated with rotated xnm
                locnew.ynm = ynmrotAll; % updated with rotated ynm
                locnew.znm = znmAll;    % updated with rotated znm
                locnew.channel = channelAll;    % original channels
                locnew.class = classAll; % original ID of ROIs
                locnew.origin = allIndLoc; % original locs
                locnew.order = orderAll; % original order of ROIs
                locnew.filenumber = repelem(newfile, length(allIndLoc))'; % assign a new file number

                 locc=obj.locData.copy;
                 locc.loc=locnew;
                 locc.regroup;
                 locc.filter;

                if p.addfile
                    locc.loc.xnm=locc.loc.xnm+x0;
                    locc.loc.ynm=locc.loc.ynm+y0;
                    fn=fieldnames(locc.loc);
                    for k=length(fn):-1:1
                        obj.locData.addloc(fn{k},locnew.(fn{k}))
                    end
                    obj.locData.regroup;
                    obj.locData.filter;
                end
           
            end
            obj.locData.regroup;
            obj.locData.filter;
            
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