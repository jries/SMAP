classdef Fisheyeview<interfaces.DialogProcessor&interfaces.SEProcessor
% Positions the sites in a fish eye view with central sites magnified
    methods
        function obj=Fisheyeview(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'se_viewer','sr_pos','sr_layeron','sr_size','se_siteroi'};
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
            x0=p.sr_pos(1);
            y0=p.sr_pos(2);
            if ~isfield(locnew,'class')
                locnew.class=0*locnew.xnm;
                locnew.filenumbernew=0*locnew.xnm;
            end
            ticc=tic;
            for k=1:length(sites)
                if  sites(k).annotation.use
                    [locs,indloc]=obj.locData.getloc({'xnm','ynm'},'position',sites(k),'grouping','ungrouped','layer',find(p.sr_layeron));
                    [phisite,rsite]=cart2pol(sites(k).pos(1)-p.sr_pos(1),sites(k).pos(2)-p.sr_pos(2));
                    maxr=mean(p.sr_size)/2; %to scale r between 0 and 1
%                     rsites=(rsite+p.offset)/maxr;
%                     rref0=p.se_siteroi(1)/maxr/12;
%                     rfac=max(rref0,rsites-(p.se_siteroi(1)/3)/maxr);
%                     rsiten=rsites.^p.scaleexponent-(p.offset/maxr).^p.scaleexponent;
                    rsiten=log(rsite+p.offset)-log(p.offset);
                    dr=1/(rsite+p.offset);
%                     maxr=1;
%                     dr=p.scaleexponent.*rfac.^(p.scaleexponent-1);
%                     dr0=p.scaleexponent.*rref0.^(p.scaleexponent-1);
%                     factor=min(dr0,dr);
                    factor=dr*1.4;
                    
                    facexphot=1.5;
                    [xsnew,ysnew]=pol2cart(phisite,rsiten);
                    locnew.xnm(indloc)=((locs.xnm-sites(k).pos(1))*factor+xsnew)*maxr;
                    locnew.ynm(indloc)=((locs.ynm-sites(k).pos(2))*factor+ysnew)*maxr;
                    locnew.phot(indloc)=locnew.phot(indloc)*factor^facexphot;
                    locnew.locprecnm(indloc)=max(10,locnew.locprecnm(indloc))*factor*maxr;
                    %also locprecnm
                    %rescale photons by factor^2
    %                 figure(88)
    %                 plot(locnew.xnm(indloc),locnew.ynm(indloc),'+')
                    locnew.filenumber(indloc)=newfile;
                    locnew.class(indloc)=sites(k).ID;
                    locnew.filenumbernew(indloc)=newfile;
                    used=used|indloc;
                    drall(k)=dr;
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
%              locc.regroup;
             locc.filter;

            if p.addfile
                locnew.xnm=locnew.xnm+x0;
                locnew.ynm=locnew.ynm+y0;
                 obj.locData.addLocData(locnew);
%                 for k=1:length(fn)
%                     obj.locData.addloc(fn{k},locnew.(fn{k}))
%                 end
                %directly modify grouploc. regroup doesnt work because of
                %change
%                 obj.locData.regroup;
                obj.locData.filter;
            end
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end




function pard=guidef

pard.t1.object=struct('String','Fisheye view. Center from main SMAP reconstruction window','Style','text');
pard.t1.position=[1,1];
pard.t1.Width=4;

pard.scaleexponentt.object=struct('String','Scaling exponent','Style','text');
pard.scaleexponentt.position=[2,1];
pard.scaleexponentt.Width=1.5;

pard.scaleexponent.object=struct('String','0.4','Style','edit');
pard.scaleexponent.position=[2,2.5];
pard.scaleexponent.Width=.5;

pard.offsett.object=struct('String','offset R (nm)','Style','text');
pard.offsett.position=[2,3];
pard.offsett.Width=1.5;

pard.offset.object=struct('String','0','Style','edit');
pard.offset.position=[2,4.5];
pard.offset.Width=.5;

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
pard.plugininfo.description='Positions the sites in a fish eye view with central sites magnified';

end