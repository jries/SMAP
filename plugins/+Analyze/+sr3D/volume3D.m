classdef volume3D<interfaces.DialogProcessor
    % volume3D renders dataset as 3D volumes with 3D elliptical Gaussiancs
    % corresponding to locprecnm and locprecznm.
    properties
        imagestack
        imageinfo
    end
    methods
        function obj=volume3D(varargin)        
            obj@interfaces.DialogProcessor(varargin{:}) ;
            obj.inputParameters={'sr_pixrec','linewidth_roi','layer1_'};
            obj.showresults=true;
        end
        function makeGui(obj)
            makeGui@interfaces.DialogProcessor(obj);
        end
        function out=run(obj,p)
            make3Dvolume(obj,p)           
            out=[];
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function make3Dvolume(obj,p)
global SMAP_stopnow
lochere=obj.locData.copy;
[locsout,indout,hroi]=lochere.getloc({'xnm','ynm','znm','xnmline','ynmline'},'position','roi');
lochere.removelocs(~indout);

if ~isempty(locsout.xnmline)
    lochere.loc.xnm=locsout.xnmline;
    lochere.loc.ynm=locsout.ynmline;
    phere.sr_pos=[0 0];
    pos=hroi.getPosition;
    dx=(pos(2,1)-pos(1,1))*1000;
    dy=(pos(2,2)-pos(1,2))*1000;
    len=sqrt(dx^2+dy^2);
   
    phere.sr_size=[len p.linewidth_roi]/2;
elseif isvalid(hroi)
    pos=hroi.getPosition;
    phere.sr_pos=[pos(1)+pos(3)/2,pos(2)+pos(4)/2]*1000;
    phere.sr_size=[pos(3) pos(4)]*1000/2;
end

if p.pixzauto
    pixz=p.pixzrecset;
else
    pixz=p.sr_pixrec;
end

if p.pixxyauto
    pixxy=p.pixxyrecset;
else
    pixxy=p.sr_pixrec;
end

obj.imageinfo.pixelsizex=pixxy;
obj.imageinfo.pixelsizez=pixz;
lochere.regroup;
lochere.filter;

 pall=obj.getLayerParameters;

 phere.sr_axes=[];
 phere.sr_pixrec=pixxy;
 
 
for k=1:length(pall)
    pall{k}=copyfields(pall{k},p);
    pall{k}=copyfields(pall{k},phere);
end
 phere.imaxtoggle=1;
  phere.imax_min=1;
imall=TotalRender(lochere,pall,{'xnm','ynm'});
 phere.imaxtoggle=0;
 phere.imax_min=imall.imax;
 
 for k=1:length(pall)
    pall{k}=copyfields(pall{k},p);
    pall{k}=copyfields(pall{k},phere);
end
sim=size(imall.image);
ax=obj.initaxis('image');

z=p.zmin:pixz:p.zmax;
imout3=zeros(sim(1),sim(2),3,length(z));

showtic=tic;
updatetime=1;
% if p.sigmazauto
%     sz=lochere.loc.znm*0+p.sigmaz;
% else
    sz=max((lochere.loc.locprecznm)*p.layer1_.gaussfac,p.layer1_.mingausspix*pixz);%gaussfac etc
    szg=max((lochere.grouploc.locprecznm)*p.layer1_.gaussfac,p.layer1_.mingausspix*pixz);
% end

for k=1:length(z)
    dz=lochere.loc.znm-z(k);
    dzg=lochere.grouploc.znm-z(k);
    lochere.loc.intensity_render=exp(-(dz.^2/2./sz.^2))/sqrt(2*pi)./sz*pixz;
  lochere.grouploc.intensity_render=exp(-(dzg.^2/2./szg.^2))/sqrt(2*pi)./szg*pixz;

    srim=TotalRender(lochere,pall,{'xnm','ynm'});
    imout3(:,:,:,k)=imout3(:,:,:,k)+srim.composite;
    imdraw=srim.composite;
    imdraw=imdraw/max(imdraw(:));
    if toc(showtic)>updatetime
        imagesc(imdraw,'Parent',ax);
        title(z(k))
        axis(ax,'equal')
         drawnow
         showtic=tic;
    end
    if SMAP_stopnow
        break
    end
end
imagesc(imall.composite/max(imall.composite(:)),'Parent',ax);
improject=permute(squeeze(sum(imout3,1)),[3 1 2]);

ax2=obj.initaxis('projection');
imagesc(improject/max(improject(:)),'Parent',ax2);
axis(ax2,'equal')
obj.imagestack=imout3;
end

function save_callback(a,b,obj)
options.color=true;
options.message=true;
options.comp='lzw';
[p,f]=fileparts(obj.locData.files.file(1).name);
[f,p]=uiputfile([p filesep f '_3Dvolume.tif']);
fileout=[p f];
if f
    imout3=obj.imagestack;
    imout=uint8(imout3/max(imout3(:))*(2^8-1));
    saveastiff(imout,fileout,options)
end
end

function openfiji_callback(a,b,obj)
title=obj.getPar('layer1_').ch_filelist.selection;
openstackinfiji(obj,obj.imagestack,title)
end

function volumeviewer_callback(a,b,obj)
imV=squeeze(sum(obj.imagestack,3));
p=obj.getGuiParameters;
imV=imV-min(imV(:));

if p.filter3dc
    
  sx=p.filter3dv(1)/obj.imageinfo.pixelsizex;
  ssx=round(5*sx);
  h=fspecial('gauss',ssx,sx);
  b=h(round(end/2),:);
  b=b/sum(b);
  
  imVx=filter(b,1,imV,[],1);
    imVy=filter(b,1,imVx,[],2);
    
  sz=p.filter3dv(end)/obj.imageinfo.pixelsizez;
  ssz=round(3*sz);
  h=fspecial('gauss',ssz,sz);
  b=h(round(end/2),:);
   b=b/sum(b);
  imV=filter(b,1,imVy,[],3);
 
end
 imV=imV/max(imV(:));
% if p.invert
%     imV=1-imV;
% end
volumeViewer(imV);
end

function pard=guidef(obj)
pard.text2.object=struct('String','zmin','Style','text');
pard.text2.position=[2,1];
pard.text3.object=struct('String','zmax','Style','text');
pard.text3.position=[3,1];

pard.zmin.object=struct('Style','edit','String','-400'); 
pard.zmin.position=[2,2.5];
pard.zmax.object=struct('Style','edit','String','400'); 
pard.zmax.position=[3,2.5];

pard.pixxyauto.object=struct('Style','checkbox','String','set pixelsize in xy (nm) to:','Value',0);
pard.pixxyauto.position=[4,1];
pard.pixxyauto.Width=1.5;
pard.pixxyrecset.object=struct('Style','edit','String','5'); 
pard.pixxyrecset.position=[4,2.5];

pard.pixzauto.object=struct('Style','checkbox','String','set pixelsize in z (nm) to:','Value',0);
pard.pixzauto.position=[5,1];
pard.pixzauto.Width=1.5;
pard.pixzrecset.object=struct('Style','edit','String','5'); 
pard.pixzrecset.position=[5,2.5];

pard.sigmazauto.object=struct('String','set sigma_z (nm) to:','Style','checkbox');
pard.sigmazauto.position=[6,1];
pard.sigmazauto.Width=1.5;
pard.sigmaz.object=struct('Style','edit','String','20'); 
pard.sigmaz.position=[6,2.5];


pard.filter3dc.object=struct('String','3D filter sIgma (nm):','Style','checkbox');
pard.filter3dc.position=[7,2.5];
pard.filter3dc.Width=1.5;
pard.filter3dv.object=struct('Style','edit','String','15'); 
pard.filter3dv.position=[7,4];
% pard.invert.object=struct('Style','checkbox','String','invert'); 
% pard.invert.position=[3,3];

pard.save.object=struct('Style','pushbutton','String','Save','Callback',{{@save_callback,obj}});
pard.save.position=[2,4];

pard.openfiji.object=struct('Style','pushbutton','String','open in Fiji','Callback',{{@openfiji_callback,obj}});
pard.openfiji.position=[3,4];

pard.volumeviewer.object=struct('Style','pushbutton','String','Matlab volumeViewer','Callback',{{@volumeviewer_callback,obj}});
pard.volumeviewer.position=[7,1];
pard.volumeviewer.Width=1.5;

pard.plugininfo.name='3D volume';
pard.plugininfo.type='ProcessorPlugin';

pard.plugininfo.description= 'volume3D renders dataset as 3D volumes with 3D elliptical Gaussiancs corresponding to locprecnm and locprecznm.';
end
