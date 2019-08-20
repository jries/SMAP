classdef sideview<interfaces.DialogProcessor
    % reconstructs side view from localizations in selected ROI
    methods
        function obj=sideview(varargin)        
            obj@interfaces.DialogProcessor(varargin{:}) ;
            obj.inputParameters=anyRender;
            obj.inputParameters{end+1}='sr_roihandle';
            obj.inputParameters{end+1}='linewidth_roi';
             obj.showresults=true;

        end
        
        function out=run(obj,p)
            sideview_reconstruct(obj,obj.locData.copy,p);  
            out=0;
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end

function sideview_reconstruct(obj,locData,p)
[locs,indroi]=locData.getloc({'xnm','ynm','locprecnm','xnmline','ynmline'},'grouping','ungrouped','position','roi');
[glocs,gindroi]=locData.getloc({'xnm','ynm','locprecnm','xnmline','ynmline'},'grouping','grouped','position','roi');

 locData.loc.xnmp=0*locData.loc.xnm-inf;
 locData.grouploc.xnmp=0*locData.grouploc.xnm-inf;
 locData.loc.ynmp=0*locData.loc.ynm-inf;
 locData.grouploc.ynmp=0*locData.grouploc.ynm-inf;
 
if ~isempty(locs.xnmline) %isfield(locs,'xnmline')
    locData.loc.xnmp(indroi)=locs.xnmline;
    locData.grouploc.xnmp(gindroi)=glocs.xnmline;
    roih=p.sr_roihandle;
    if isvalid(roih)
        rpos=roih.getPosition;
        lr=sqrt(sum((rpos(2,:)-rpos(1,:)).^2));
        rx=[-lr/2 lr/2]*1000;
        ry=[-0.5 0.5]*p.linewidth_roi;
    else
        rx=[min(locs.xnmline) max(locs.xnmline)];
        ry=[min(locs.ynmline) max(locs.ynmline)];
    end
    locData.loc.ynmp(indroi)=locs.ynmline;
    locData.grouploc.ynmp(gindroi)=glocs.ynmline;
else
    locData.loc.xnmp(indroi)=locs.xnm;
    locData.grouploc.xnmp(gindroi)=glocs.xnm; 
    locData.loc.ynmp(indroi)=locs.ynm;
    locData.grouploc.ynmp(gindroi)=glocs.ynm; 
    rx=[min(locs.xnmp) max(locs.xnmp)];
    ry=[min(locs.ynmp) max(locs.ynmp)];
end

if ~p.pixauto
    pixelsizex=p.sr_pixrec;pixelsizey=p.sr_pixrec;pixelsizez=p.sr_pixrec;
else
    pixelsizex=p.pixrecset(1);pixelsizey=p.pixrecset(1);pixelsizez=p.pixrecset(end);
end
rz=[p.zmin p.zmax];

yax=obj.initaxis('x-y');
p.sr_axes=[];p.rangex=rx;p.rangey=ry;p.sr_pixrec=[pixelsizex pixelsizey];
if isfield(locData.loc,'locprecznm')
    fieldsz='locprecznm';
else
    fieldsz='locprecnm';
end
yim=anyRender(locData,p,'x','xnmp','y','ynmp','sx','locprecnm','sy',fieldsz);
imagesc(rx,ry,yim.image,'Parent',yax);
axes(yax)
axis equal;
axis tight
p.sr_axes=[];
zax=obj.initaxis('x-z');
p.rangex=rx;p.rangey=rz;p.sr_pixrec=[pixelsizex pixelsizez];
[zim]=anyRender(locData,p,'x','xnmp','y','znm','sx','locprecnm','sy',fieldsz);
imagesc(rx,rz,zim.image,'Parent',zax);
% [~,ax1]=renderAnyField(locData,p,'x-z');
axes(zax)
axis equal;
axis ij;
axis tight;
s=size(yim.image);
border=ones(1,s(2),3);
comp=vertcat(yim.image,border,zim.image);
rz(1)=rz(1)-(ry(2)-ry(1));
cax=obj.initaxis('x-z/x-y');
imagesc(rx,rz,comp,'Parent',cax);
axis equal;
axis ij;
axis tight
end


function pard=guidef
pard.text1.object=struct('String','parameters','Style','text');
pard.text1.position=[1,1];

pard.text2.object=struct('String','zmin','Style','text');
pard.text2.position=[2,1];
pard.text3.object=struct('String','zmax','Style','text');
pard.text3.position=[3,1];

pard.zmin.object=struct('Style','edit','String',-400); 
pard.zmin.position=[2,2];
pard.zmax.object=struct('Style','edit','String',400); 
pard.zmax.position=[3,2];

pard.pixauto.object=struct('Style','checkbox','String','set pixelsize (x,z)','Value',0);
pard.pixauto.position=[4,1];
pard.pixrecset.object=struct('Style','edit','String','5, 5'); 
pard.pixrecset.position=[4,2];

pard.plugininfo.name='Side view';
pard.plugininfo.description= 'Reconstructs side view from localizations in selected ROI';
pard.plugininfo.type='ProcessorPlugin';
end