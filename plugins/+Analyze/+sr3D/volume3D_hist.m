classdef volume3D_hist<interfaces.DialogProcessor
    % volume3D renders dataset as 3D volumes with 3D elliptical Gaussiancs
    % corresponding to locprecnm and locprecznm.
    properties
        imagestack
        imageinfo
        fileinfo
    end
    methods
        function obj=volume3D_hist(varargin)        
            obj@interfaces.DialogProcessor(varargin{:}) ;
            obj.inputParameters={'sr_pixrec','linewidth_roi','layer1_','sr_layerson'};
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

phere=p;

layers=find(p.sr_layerson);
lochere=obj.locData;
for k=1:length(layers) 
    [locsout{k},indout,hroi]=lochere.getloc({'xnm','ynm','znm','xnmline','ynmline'},'position','roi','layer',layers(k));
end


if ~isempty(locsout{1}.xnmline)
    for k=1:length(layers)
        locsout{k}.xnm=locsout{k}.xnmline;
        locsout{k}.ynm=locsout{k}.ynmline;
    end
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

minx=phere.sr_pos(1)-phere.sr_size(1);maxx=phere.sr_pos(1)+phere.sr_size(1);
miny=phere.sr_pos(2)-phere.sr_size(2);maxy=phere.sr_pos(2)+phere.sr_size(2);
minz=phere.zmin;maxz=phere.zmax;

if p.uint8
    format='uint8';
else
    format='uint16';
end
   

sizeV=ceil([(maxx-minx)/pixxy (maxy-miny)/pixxy (maxz-minz)/pixz]);
imout=zeros(sizeV(1),sizeV(2),sizeV(3),length(layers),format);

for k=1:length(layers)
    xh=ceil((locsout{k}.xnm-minx)/pixxy);
    yh=ceil((locsout{k}.ynm-miny)/pixxy);
    zh=ceil((locsout{k}.znm-minz)/pixz);
    
    indg=xh>0&xh<=sizeV(1) & yh>0&yh<=sizeV(2) & zh>0&zh<=sizeV(3);
    
    linind=sub2ind(sizeV,xh(indg),yh(indg),zh(indg));
    hl=cast(histc(linind,1:max(linind)),format);
%     hl=cast(histcounts(linind,1:max(linind)),format);
    disp(['maxcounts:' num2str(max(hl))])
    uind=unique(linind);
    imouth=zeros(sizeV,format);
    imouth(uind)=hl(hl>0);
    imout(:,:,:,k)=imouth;
end

%save txtfile with info:
obj.fileinfo.position=pos;
obj.fileinfo.type=class(hroi);
obj.fileinfo.lienwidth=p.linewidth_roi;
obj.fileinfo.minz=minz;
obj.fileinfo.maxz=maxz;
obj.fileinfo.sr_pos=phere.sr_pos;
obj.fileinfo.sr_size=phere.sr_size;
obj.fileinfo.pixxy=pixxy;
obj.fileinfo.pixz=pixz;
% imageslicer(imout);
obj.imagestack=imout;
end

function save_callback(a,b,obj)
options.color=false;
options.message=true;
options.comp='lzw';
[p,f]=fileparts(obj.locData.files.file(1).name);
[f,p]=uiputfile([p filesep f '_3Dvolume.tif']);
fileout=[p f];
if f
    imout=obj.imagestack;
    for k=1:size(imout,4)
        fileouth=strrep(fileout,'.tif',['_' num2str(k) '.tif']);
%     imout=uint8(imout3/max(imout3(:))*(2^8-1));
    saveastiff(imout(:,:,:,k),fileouth,options)
    end
    
    fileoutt=strrep(f,'.tif','.txt');
    writestruct([p fileoutt],obj.fileinfo);
    
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


pard.uint8.object=struct('Style','checkbox','String','use 8 bit'); 
pard.uint8.position=[6,1];

% pard.sigmazauto.object=struct('String','set sigma_z (nm) to:','Style','checkbox');
% pard.sigmazauto.position=[6,1];
% pard.sigmazauto.Width=1.5;
% pard.sigmaz.object=struct('Style','edit','String','20'); 
% pard.sigmaz.position=[6,2.5];


% pard.filter3dc.object=struct('String','3D filter sIgma (nm):','Style','checkbox');
% pard.filter3dc.position=[7,2.5];
% pard.filter3dc.Width=1.5;
% pard.filter3dv.object=struct('Style','edit','String','15'); 
% pard.filter3dv.position=[7,4];
% pard.invert.object=struct('Style','checkbox','String','invert'); 
% pard.invert.position=[3,3];

pard.save.object=struct('Style','pushbutton','String','Save','Callback',{{@save_callback,obj}});
pard.save.position=[2,4];

pard.openfiji.object=struct('Style','pushbutton','String','open in Fiji','Callback',{{@openfiji_callback,obj}});
pard.openfiji.position=[3,4];

pard.volumeviewer.object=struct('Style','pushbutton','String','Matlab volumeViewer','Callback',{{@volumeviewer_callback,obj}});
pard.volumeviewer.position=[7,1];
pard.volumeviewer.Width=1.5;

pard.plugininfo.name='3D_volume_hist';
pard.plugininfo.type='ProcessorPlugin';

pard.plugininfo.description= 'volume3D renders dataset as 3D volumes with 3D elliptical Gaussiancs corresponding to locprecnm and locprecznm.';
end
