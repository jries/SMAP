classdef segmentBeads<interfaces.DialogProcessor&interfaces.SEProcessor
%     Segements beads
    methods
        function obj=segmentBeads(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'se_cellfov','se_sitefov','se_siteroi'};
        end
        
        function out=run(obj,p)  
          segmentb(obj,p);
          out=[];
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end


function segmentb(obj,p)

locDatacopy=obj.locData.copy;
locDatacopy.regroup(150,p.minframes);
% g=Grouper(locDatacopy);
% g.connect(150,max(locDatacopy.loc.frame),'frame','xnm','ynm','filenumber');
% g.connect(50,3,'frame','xnm','ynm','filenumber');
% g.combine;
beadind=locDatacopy.grouploc.numberInGroup>p.minframes;
xb=locDatacopy.grouploc.xnm(beadind);
yb=locDatacopy.grouploc.ynm(beadind);
filenumber=locDatacopy.grouploc.filenumber(beadind);
% p=obj.getAllParameters;
% files=obj.SE.files;
% xm=min(obj.locData.loc.xnm);
% ym=min(obj.locData.loc.ynm);
% xx=max(obj.locData.loc.xnm);
% yx=max(obj.locData.loc.ynm);
% posx=xm-p.se_sitefov:p.se_cellfov:xx+p.se_sitefov;
% posy=ym-p.se_sitefov:p.se_cellfov:yx+p.se_sitefov;
% for k=1:length(files)
%     for x=1:length(posx)
%         for y=1:length(posy)
%             pos=[posx(x) posy(y)]+p.se_cellfov/2;
% 
%              currentcell=interfaces.SEsites;
%             currentcell.pos=pos;
%             currentcell.ID=0;
% 
%             currentcell.info.filenumber=k;
%             obj.SE.addCell(currentcell);
cells=obj.SE.cells;
for c=1:length(cells)
    cellh=cells(c);
    pos=cellh.pos;
    cfn=cellh.info.filenumber;
            beadhere=(mywithin(xb,[pos(1)-p.se_cellfov/2 p.se_cellfov],yb,[pos(2)-p.se_cellfov/2 p.se_cellfov]));
            beadfile=filenumber==cfn;
            beadhere=beadhere&beadfile;
            beadhere=find(beadhere);
            for b=1:length(beadhere)
                    currentsite=interfaces.SEsites;
                    currentsite.pos=[xb(beadhere(b)) yb(beadhere(b)) 0];     
                    currentsite.info.cell=cellh.ID;
                    currentsite.info.filenumber=cellh.info.filenumber;
                     obj.SE.addSite(currentsite);
            end
%         end
%     end
end
obj.SE.processors.preview.updateCelllist
obj.SE.processors.preview.updateSitelist
end

function filterbeads(a,b,obj)
p=obj.getAllParameters;
ax=obj.initaxis('sx');
hold on
ax2=obj.initaxis('sy');
hold on
axp=obj.initaxis('phot');
hold on
sites=obj.SE.sites;
framezero=(max(obj.locData.loc.frame)+min(obj.locData.loc.frame))/2;
framerange=10;
for s=length(sites):-1:1
    if sites(s).annotation.use
        posh=horzcat(sites(s).pos(1:2),p.se_siteroi);
        locs=obj.locData.getloc({'xnm','ynm','PSFxnm','PSFynm','frame','phot'},'Position',posh,'layer',1,'grouping','ungrouped');
        
        plot(ax, locs.frame,locs.PSFxnm,'y-')
        hold on
         plot(ax2, locs.frame,locs.PSFynm,'y-')
        hold on       
        plot(axp, locs.frame,locs.phot,'.')
        hold on
        
        if length(p.framerange)>1
            range=locs.frame>p.framerange(1)&&locs.frame<p.framerange(2);
        else
        range=find(abs(locs.frame-framezero)<p.framerange);
        end
        psfx=locs.PSFxnm;
        psfy=locs.PSFynm;
        [minPSFx(s),indx]=min(psfx(range));
        [minPSFy(s),indy]=min(psfy(range));
        
%         [~,indx]=min(locs.PSFxnm./(locs.phot-min(locs.phot)));
%         [~,indy]=min(locs.PSFynm./(locs.phot-min(locs.phot)));
        
        frameminx(s)=locs.frame(range(indx));
        frameminy(s)=locs.frame(range(indy));
        [maxphot(s),indp]=max(locs.phot);
        framephot(s)=locs.frame(indp);
    end
    
end
if ~isempty(obj.SE.currentsite)
posh=horzcat(obj.SE.currentsite.pos(1:2),p.se_siteroi);
locs=obj.locData.getloc({'xnm','ynm','PSFxnm','PSFynm','frame','phot'},'Position',posh,'layer',1,'grouping','ungrouped');

plot(ax, locs.frame,locs.PSFxnm,'k-')
hold on
plot(ax2, locs.frame,locs.PSFynm,'k-')
hold on
plot(axp, locs.frame,locs.phot,'k-')
hold on
end
axh=obj.initaxis('PSFxpos')
plot(frameminx,minPSFx,'.')
hold on
% axh=obj.initaxis('histmin')
plot(frameminy,minPSFy,'.')
hold off
axh=obj.initaxis('photf')
plot(framephot,maxphot,'.')
psfdist=p.psfdist;
mindist=p.framedistance;
if p.psfdistc
indgood=(minPSFx)<psfdist+median(minPSFx(minPSFx>0))&(minPSFy)<psfdist+median(minPSFy(minPSFy>0));
else
    indgood=true(size(minPSFx));
end
if p.framedistancec
 indgood=indgood&abs(frameminx-median(frameminx(frameminx>0)))<=mindist&abs(frameminy-median(frameminy(frameminy>0)))<=mindist;
end
remove=find(~indgood);
for s=1:length(remove)
    
    sites(remove(s)).annotation.use=false;
end
obj.SE.processors.preview.updateSitelist


for s=length(sites):-1:1
    if sites(s).annotation.use
        posh=horzcat(sites(s).pos(1:2),p.se_siteroi);
        locs=obj.locData.getloc({'xnm','ynm','PSFxnm','PSFynm','frame','phot'},'Position',posh,'layer',1,'grouping','ungrouped');
        
        plot(ax, locs.frame,locs.PSFxnm,'-')
        hold on
        plot(ax2, locs.frame,locs.PSFynm,'-')
        hold on        
        plot(axp, locs.frame,locs.phot,'.')
        hold on
       
    end
    
end

end

function makecells(a,b,obj)
f=figure;%f.Visible='off';
mc=plugin('ROIManager','Segment','makeCellGrid',f,obj.P);
mc.attachLocData(obj.locData)
mc.attachSE(obj.locData.SE)
mc.makeGui;
p=mc.getAllParameters;
mc.run(p);
close(f)
end

function pard=guidef(obj)




pard.makecells.object=struct('String','Make cell grid','Style','pushbutton','Callback',{{@makecells, obj}});
pard.makecells.position=[1,1];


% pard.segmentbeads.object=struct('String','Make cell grid','Style','pushbutton','Callback',{{@makecells, obj}});
% pard.segmentbeads.position=[2,1];

pard.t2.object=struct('String','minimum frames','Style','text');
pard.t2.position=[2,1];
pard.minframes.object=struct('String','30','Style','edit');
pard.minframes.position=[2,2];

pard.dofilter.object=struct('String','Filter','Style','pushbutton','Callback',{{@filterbeads, obj}});
pard.dofilter.position=[3,1];

pard.psfdistc.object=struct('String','max PSF distance from median (nm)','Style','checkbox');
pard.psfdistc.position=[4,1];
pard.psfdistc.Width=2;
pard.psfdist.object=struct('String','10','Style','edit');
pard.psfdist.position=[4,3];

pard.framedistancec.object=struct('String','max PSF distance from median (frames)','Style','checkbox');
pard.framedistancec.position=[5,1];
pard.framedistancec.Width=2;
pard.framedistance.object=struct('String','2','Style','edit');
pard.framedistance.position=[5,3];

pard.frameranget.object=struct('String','frame range considered','Style','text');
pard.frameranget.position=[6,1];
pard.frameranget.Width=2;
pard.framerange.object=struct('String','15','Style','edit');
pard.framerange.position=[6,3];
pard.plugininfo.description='Segements beads';
% pard.t3.object=struct('String','diameterNPC','Style','text');
% pard.t3.position=[3,1];
% pard.diameterNPC.object=struct('String','110','Style','edit');
% pard.diameterNPC.position=[3,2];
% pard.t4.object=struct('String','rim','Style','text');
% pard.t4.position=[4,1];
% pard.rim.object=struct('String','20','Style','edit');
% pard.rim.position=[4,2];
% 
% pard.saveon.object=struct('String','saveon','Style','checkbox');
% pard.saveon.position=[1,3];
% 
% pard.getmask.object=struct('String','getmask','Style','checkbox');
% pard.getmask.position=[2,3];
% pard.plugininfo.type='ROI_Analyze';
end

