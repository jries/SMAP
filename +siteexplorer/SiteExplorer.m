classdef SiteExplorer<rec.LocProcessor
    properties
        sePar
        sites
        numberOfSites=0;
        maxsite=0;
        cells
        numberOfCells=0;
        numberOfFiles=0;
        maxcell=0;
        processors
        files
        currentsiteID=0;
        currentcellID=0;
        currentfileID=1;
        currentsite
        currentcell
        currentfile
        
        temp
%         sitetemp
%         celltemp
%         filetemp
    end
    methods
        function obj=SiteExplorer(varargin)
            obj.sites=siteexplorer.SEsites();
            obj.cells=siteexplorer.SEsites();
            obj.files=siteexplorer.SEsites();
            obj.temp.siteincell=[];
            obj.temp.cellinfile=[];
%             obj.cellboxes=siteexplorer.boximage;
%             obj.siteboxes=siteexplorer.boximage;

        end
        function empty(obj)
            obj.sites=siteexplorer.SEsites();
            obj.cells=siteexplorer.SEsites();
            obj.files=siteexplorer.SEsites();
            obj.numberOfSites=0;
            obj.maxsite=0;
            obj.numberOfCells=0;
            obj.numberOfFiles=0;
            obj.maxcell=0;
            obj.currentsiteID=0;
            obj.currentcellID=0;
        end
        function attachLocData(obj,locData)
            attachLocData@rec.LocProcessor(obj,locData);
            obj.initializeProcessors()
        end
        function addSites(obj, SEin,filenew)
            
            %renumber files
            if ~isempty(SEin.files)
            for k=1:length(SEin.files)
                SEin.files(k).ID=filenew(SEin.files(k).ID);
            end
            for k=1:length(SEin.cells)
                SEin.cells(k).info.filenumber=filenew(SEin.cells(k).info.filenumber);
                SEin.cells(k).ID=SEin.cells(k).ID+obj.maxcell;
            end
            for k=1:length(SEin.sites)
                SEin.sites(k).info.filenumber=filenew(SEin.sites(k).info.filenumber);
                SEin.sites(k).info.cell=SEin.sites(k).info.cell+obj.maxcell;
                SEin.sites(k).ID=SEin.sites(k).ID+obj.maxsite;
            end
            
            
            %files
            SEin.numberOfFiles=length(SEin.files);
            obj.files(obj.numberOfFiles+1:obj.numberOfFiles+SEin.numberOfFiles)=SEin.files;
            obj.numberOfFiles=obj.numberOfFiles+SEin.numberOfFiles;
            
%           cells
            SEin.numberOfCells=length(SEin.cells);
            obj.cells(obj.numberOfCells+1:obj.numberOfCells+SEin.numberOfCells)=SEin.cells;
            obj.numberOfCells=obj.numberOfCells+SEin.numberOfCells;
%            sites
            SEin.numberOfSites=length(SEin.sites);
            obj.sites(obj.numberOfSites+1:obj.numberOfSites+SEin.numberOfSites)=SEin.sites;
            obj.numberOfSites=obj.numberOfSites+SEin.numberOfSites;            
            
            obj.maxsite=max([obj.sites.ID]);
            obj.maxcell=max([obj.cells.ID]);
            
            
            obj.currentfile=obj.files(1);
            obj.currentcell=obj.cells(1);
            obj.currentsite=obj.sites(1);
            obj.setIndList;
            end
            
        end
        function setIndList(obj)
            for k=1:length(obj.sites)
                obj.sites(k).indList=k;
            end
        end
        function ID=addSite(obj,site)
            if (isempty(site.ID)||site.ID==0)
                obj.numberOfSites=obj.numberOfSites+1;
                obj.maxsite=obj.maxsite+1;
                site.ID=obj.maxsite;
                site.indList=obj.numberOfSites;
                obj.sites(obj.numberOfSites)=site;
                
            else 
                disp('site already there, ID not 0')
                
            end
            ID=site.ID;
        end
        function removeSite(obj,siteID)
            ind=obj.indexFromID(obj.sites,siteID);
            obj.numberOfSites=obj.numberOfSites-1;
            obj.sites(ind)=[];
        end
        function removeCell(obj,cellID)
            ind=obj.indexFromID(obj.cells,cellID);
            obj.cells(ind)=[];
            si=[obj.sites.info];
            indout=[si.cell]==cellID;
            obj.sites(indout)=[];
            obj.numberOfCells=obj.numberOfCells-1;
            obj.numberOfSites=obj.numberOfSites-sum(indout);
        end
        function ID=addCell(obj,cell)
            obj.numberOfCells=obj.numberOfCells+1;
            obj.maxcell=obj.maxcell+1;
            cell.ID=obj.maxcell;
            obj.cells(obj.numberOfCells)=cell;
            ID=cell.ID;
            
%             %boxes
%             cellsize=obj.sePar.Settings.cellfov;
%             if isempty(obj.cellboxes.boxsizenm)||obj.cellboxes.boxsizenm~=cellsize
%                 obj.cellboxes.init(
%             end
            
            
        end
        function addFile(obj,name,ID,info)
            ind=obj.indexFromID(obj.files,ID);
            if isempty(ind)
                obj.numberOfFiles=obj.numberOfFiles+1;
                ind=obj.numberOfFiles;
            end
            
            file=siteexplorer.SEsites;
            file.name=name;
            file.ID=ID;
            file.info=info;
            obj.files(ind)=file;
        end
        
        function initializeProcessors(obj)
            obj.processors.renderer=rec.Renderer(obj.locData);
            obj.processors.drawer=rec.Drawer(obj.locData);
            obj.processors.displayer=rec.Displayer(obj.locData);
        end
        
        function imout=plotsite(obj,site,hax,hbox)
            
            if isnumeric(site)
                ind=obj.indexFromID(obj.sites,site);
                site=obj.sites(ind);
            end
            
            if isempty(site.image)
                
%                  display('draw site')
%             p1=obj.locData.parameters;
            p1.sr_pos=site.pos;
            p1.sr_size=ones(2,1)*obj.sePar.Settings.sitefov/2;
            p1.pixrec=obj.sePar.Settings.sitepixelsize;
            p1.sr_sizeRecPix=round((p1.sr_size*2)/p1.pixrec);
            p1.sr_axes=hax;
            p1.sr_axes=-1;
            p1.normalizeFoV=p1.sr_sizeRecPix(1)/obj.sePar.Settings.sitefov*obj.sePar.Settings.siteroi/2;
            angle=siteexplorer.pos2angle(site.annotation.rotationpos.pos);
            if obj.sePar.Settings.rotate&&angle~=0
             p1.rotationangle=angle;
            else 
             p1.rotationangle=0;
            end
            site.image=obj.plotobject(p1,site.info.filenumber);%filenumber
            site.image.angle=p1.rotationangle; %remove later? not needed
%             site.handles.roi=hggroup('Parent',hax);
% %             plotbox(site.handles.roi,site.pos,obj.sePar.Settings.siteroi)
%             plotbox(hax,site.pos,obj.sePar.Settings.siteroi)
            end
            
             displayimage(site.image,hax);
             plotbox(hax,site.pos,obj.sePar.Settings.siteroi);
             plotcirc(hax,site.pos,obj.sePar.Settings.siteroi);
             delete(obj.temp.siteincell);
             obj.temp.siteincell=plotbox(hbox,site.pos,obj.sePar.Settings.sitefov);
             
            imout=site.image;
            
        end
        
        function imout=plotcell(obj,cell,hax,hbox)
             if isempty(cell.image)
                p1.sr_pos=cell.pos;
                p1.sr_size=ones(2,1)*obj.sePar.Settings.cellfov/2;
                p1.pixrec=obj.sePar.Settings.cellpixelsize;
                p1.sr_sizeRecPix=round((p1.sr_size*2)/p1.pixrec);       
                p1.sr_axes=hax;
                p1.sr_axes=-1;
                p1.rotationangle=0;
                p1.normalizeFoV=[];
                cell.image=obj.plotobject(p1,cell.info.filenumber);%filenumber
                
             end
            displayimage(cell.image,hax)
            
            
             %plot all site boxes
             if obj.sePar.Settings.drawboxes&&~isempty(obj.sites)&&~isempty(obj.sites(1).info)
                 allsites=[obj.sites(:).info];
                 ind=[allsites.cell]==cell.ID;
                 plotmanyboxes(hax,obj.sites,ind,obj.sePar.Settings.sitefov);
             end
            
            delete(obj.temp.cellinfile);
            obj.temp.cellinfile=plotbox(hbox,cell.pos,obj.sePar.Settings.cellfov);
            imout=cell.image;
            
           
            
        end
        
        function image=plotfile(obj,fileID,hax) %same as cell: store image, temp. check
            ind=obj.indexFromID(obj.files,fileID);
            file=obj.files(ind);
            if isempty(file.image)
                pixrec=file.info.pixsize*1000;
                roi=file.info.roi*pixrec;

                % test for SML could be all wrong
                p1.sr_pos=[roi(1)+roi(3)/2 roi(2)+roi(4)/2];
                p1.sr_size=[roi(3) roi(4)]/2; 
                p1.pixrec=pixrec;
                p1.sr_sizeRecPix=round((p1.sr_size*2)/p1.pixrec);
                p1.sr_axes=-1;
                p1.mingausspix=0.8;
                p1.gaussfactor=0.1;
                p1.rotationangle=0;
                p1.normalizeFoV=[];
                file.image=obj.plotobject(p1,fileID);
                
            end
            displayimage(file.image,hax);
            
                         %plot all site boxes
             if obj.sePar.Settings.drawboxes&&~isempty(obj.cells)&&~isempty(obj.cells(1).info)
                 allcells=[obj.cells(:).info];
                 ind=[allcells.filenumber]==fileID;
                 plotmanyboxes(hax,obj.cells,ind,obj.sePar.Settings.cellfov);
             end
            
        end
        
        function image=plotobject(obj,p,filenumber)
            
           fileind=obj.indexFromID(obj.files,filenumber);
           p1=obj.locData.parameters;
           for k=1:6
                p2=obj.locData.layer(k).parameters;
                
                if p2.layercheck
                    p2.ch_filelist.value=fileind;
                    p2.ch_filelist.selection=p2.ch_filelist.string{fileind};
                p3=obj.locData.layer(k).parameters.rec_addpar;
                obj.processors.renderer.setParameters(p1,p2,p3,p)
                obj.processors.drawer.setParameters(p1,p2,p3,p);
                
                %filter filenumber
                groupc=p2.groupcheck;
                filterold=obj.locData.getFilter(k);
                filternew=filterold;
                locs=obj.locData.getloc({'filenumber','xnm','ynm'},'grouping',groupc);
                
                filternew.filenumber=(locs.filenumber==filenumber);
                filternew.xnm=rec.LocalizationFilter.minMaxFilter(locs.xnm,p.sr_pos(1)-p.sr_size(1),p.sr_pos(1)+p.sr_size(1));
                filternew.ynm=rec.LocalizationFilter.minMaxFilter(locs.ynm,p.sr_pos(2)-p.sr_size(2),p.sr_pos(2)+p.sr_size(2));
%                 filternew=myrmfield(filternew,'xnm');
%                 filternew=myrmfield(filternew,'ynm');
                obj.locData.setFilter(filternew,k)
                
                rawimage=obj.processors.renderer.render(k);
                obj.locData.setFilter(filterold,k);
                layers(k).images.finalImages=obj.processors.drawer.drawImage(rawimage);
                 layers(k).images.rawimage=rawimage;
                
           
                end
           end
           obj.processors.displayer.setParameters(p1,p2,p3,p);
           image=obj.processors.displayer.displayImage(layers);
           image.parameters=p;
           image.layers=layers;
%            axis equal
            
        end
           
        
        function ind=indexFromID(obj,list,ID)
%             ind=[];
%             if isfield(list,'ID')
            ind=find([list.ID]==ID,1,'first');
%             end
        end
        function cell=getcell(obj,ID)
            ind=obj.indexFromID(obj.cells,ID);
            cell=obj.cells(ind);
        end
        function sites=getsite(obj,ID)
            ind=obj.indexFromID(obj.sites,ID);
            sites=obj.sites(ind);
        end
        function updateSite(obj)
            if obj.currentsite.ID>0
                obj.sites(obj.indexFromID(obj.sites,obj.currentsite.ID))=obj.currentsite;
            end
        end
        function updateCell(obj)
            if obj.currentcell.ID>0
                obj.cells(obj.indexFromID(obj.cells,obj.currentcell.ID))=obj.currentcell;
            end
        end
        
        
        function se=save(obj)
            se.sePar=obj.sePar;
            se.sites=copy(obj.sites);
            for k=1:length(se.sites)
                se.sites(k).handles=[];
                se.sites(k).image=[];
            end
            se.cells=copy(obj.cells);
            for k=1:length(se.cells)
                se.cells(k).handles=[];
                se.cells(k).image=[];
            end
            
            se.files=copy(obj.files);
            for k=1:length(se.files)
                se.files(k).image=[];
                se.files(k).handles=[];
            end
%             
%             se.files.handles=[];
            se.numberOfSites=obj.numberOfSites;
            se.numberOfCells=obj.numberOfCells;
            se.numberOfFiles=obj.numberOfFiles;
        end
        function load(obj,se)
        end
    end
end

function displayimage(img,hax)
 imagesc(img.rangex,img.rangey,img.image,'Parent',hax,'Pickable','none','HitTest','off')

set(hax,'Xlim',double(img.rangex))
set(hax,'Ylim',double(img.rangey))
set(hax,'YDir','reverse')
hax.HitTest='on';
hax.PickableParts='all';
end

function hg=plotbox(h,pos,roi,color)
if nargin<4
    color=[1 1 1];
end
linewidth=1;

pos=pos/1000;
roi=roi/1000;
x1=pos(1)-roi/2;
x2=pos(1)+roi/2;
y1=pos(2)-roi/2;
y2=pos(2)+roi/2;
hg=hggroup('Parent',h);
line([x1 x1],[y1 y2],'Parent',hg,'Color',color,'LineWidth',linewidth);
line([x1 x2],[y2 y2],'Parent',hg,'Color',color,'LineWidth',linewidth);
line([x2 x2],[y2 y1],'Parent',hg,'Color',color,'LineWidth',linewidth);
line([x2 x1],[y1 y1],'Parent',hg,'Color',color,'LineWidth',linewidth);
end

function hg=plotcirc(h,pos,roi,color)
if nargin<4
    color=[1 1 1]*.5;
end
linewidth=1;

pos=pos/1000;
roi=roi/1000;
x1=pos(1)-roi/2;
x2=pos(1)+roi/2;
y1=pos(2)-roi/2;
y2=pos(2)+roi/2;
hg=rectangle('Parent',h,'Position',[x1,y1,roi,roi],'Curvature',[1,1],'EdgeColor',color,'LineStyle','--');
% hg=hggroup('Parent',h);
% line([x1 x1],[y1 y2],'Parent',hg,'Color',color,'LineWidth',linewidth);
% line([x1 x2],[y2 y2],'Parent',hg,'Color',color,'LineWidth',linewidth);
% line([x2 x2],[y2 y1],'Parent',hg,'Color',color,'LineWidth',linewidth);
% line([x2 x1],[y1 y1],'Parent',hg,'Color',color,'LineWidth',linewidth);
end


function hg=plotmanyboxes(hax,list,ind,roi)
hg=hggroup('Parent',hax);
color=[1 0 1];
for k=1:length(list)
    if ind(k)
        plotbox(hg,list(k).pos,roi,color);
    end
end
end