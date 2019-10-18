classdef MovieRunningWindow<interfaces.DialogProcessor
    % MovieRunningWindow reconstruct a movie with running average from live-cell data
    properties
        imagestack
        movie
    end
    methods
        function obj=MovieRunningWindow(varargin)        
            obj@interfaces.DialogProcessor(varargin{:}) ;
            obj.inputParameters={'sr_pixrec','linewidth_roi','layer1_','sr_layerson'};
            obj.showresults=true;
        end
        function makeGui(obj)
            makeGui@interfaces.DialogProcessor(obj);
        end
        function out=run(obj,p)
            makemovie(obj,p)           
            out=[];
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function makemovie(obj,p)
global SMAP_stopnow

lochere=obj.locData.copy;
[locsout,indout,hroi]=lochere.getloc({'xnm','ynm','znm','xnmline','ynmline','frame'},'position','roi','layer',find(p.sr_layerson),'grouping','ungrouped');
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
% elseif isvalid(hroi)
%     pos=hroi.getPosition;
%     phere.sr_pos=[pos(1)+pos(3)/2,pos(2)+pos(4)/2]*1000;
%     phere.sr_size=[pos(3) pos(4)]*1000/2;
else
    pos=hroi.getPosition;
    phere.sr_size=pos(3:4)*1000/2;
    phere.sr_pos=(pos(1:2)+pos(3:4)/2)*1000;
end

% if p.pixzauto
%     pixz=p.pixzrecset;
% else
%     pixz=p.sr_pixrec;
% end

% if p.pixxyauto
%     pixxy=p.pixxyrecset;
% else
%     pixxy=p.sr_pixrec;
% end

lochere.regroup;
lochere.filter;

 pall=obj.getLayerParameters;

 phere.sr_axes=[];
%  phere.sr_pixrec=pixxy;
 
 
for k=1:length(pall)
    pall{k}=copyfields(pall{k},p);
    pall{k}=copyfields(pall{k},phere);
end

ax=obj.initaxis('image');
% ax=gca;

minf=min(locsout.frame);maxf=max(locsout.frame);
frames=minf:p.df:maxf;
showtic=tic;
updatetime=.1;
srim=TotalRender(lochere,pall,{'frame'});
s=size(srim.composite);
imout3=zeros(s(1),s(2),s(3),length(frames)-1);

for k=length(frames)-1:-1:1
   
    minmax=[frames(k),frames(k)+p.window];
    if p.cumulative
        minmax(1)=0;
    end
    lochere.filter('frame',[],'minmax',minmax)
    srim=TotalRender(lochere,pall);
     imout3(:,:,:,k)=srim.image;
    imdraw=srim.image;
    imdraw=imdraw/max(imdraw(:));
    if toc(showtic)>updatetime
        imagesc(imdraw,'Parent',ax);
%         title(z(k))
        axis(ax,'equal')
         drawnow
         showtic=tic;
    end
     F(k)=im2frame(srim.image);
%     F(k)=getframe(ax);
    if SMAP_stopnow
        break
    end
end
obj.movie=F;
obj.imagestack=imout3;
display('done')
title('done')
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

function play_callback(a,b,obj)
persistent f
if isempty(f) || ~isvalid(f)
    f=figure;
end
figure(f)
movie(f,obj.movie,10)
end

function openfiji_callback(a,b,obj)
title=obj.getPar('layer1_').ch_filelist.selection;
openstackinfiji(obj,obj.imagestack,title)

end

function pard=guidef(obj)
pard.text2.object=struct('String','time step (frames)','Style','text');
pard.text2.position=[2,1];
pard.text3.object=struct('String','window (frames)','Style','text');
pard.text3.position=[3,1];

pard.df.object=struct('Style','edit','String','1000'); 
pard.df.position=[2,2.5];
pard.window.object=struct('Style','edit','String','5000'); 
pard.window.position=[3,2.5];

pard.cumulative.object=struct('Style','checkbox','String','cumulative'); 
pard.cumulative.position=[4,1];


pard.save.object=struct('Style','pushbutton','String','Save','Callback',{{@save_callback,obj}});
pard.save.position=[2,4];

pard.openfiji.object=struct('Style','pushbutton','String','open in Fiji','Callback',{{@openfiji_callback,obj}});
pard.openfiji.position=[3,4];
pard.play.object=struct('Style','pushbutton','String','play','Callback',{{@play_callback,obj}});
pard.play.position=[4,4];

pard.plugininfo.name='Movie';
pard.plugininfo.type='ProcessorPlugin';

pard.plugininfo.description= 'Reconstruct a movie with running average from live-cell data';
end
