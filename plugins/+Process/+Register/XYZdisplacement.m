classdef XYZdisplacement<interfaces.DialogProcessor
    properties
    end
    methods
        function obj=XYZdisplacement(varargin)        
            obj@interfaces.DialogProcessor(varargin{:});
             obj.inputParameters={'sr_pixrec','sr_pos','sr_size'};
            obj.showresults=true;
            obj.history=false;
        end
        
        function out=run(obj,p)
            out=[];
            locr=obj.locData.getloc({'xnm','ynm','znm','locprecnm','locprecznm'},'layer',p.reflayer.Value,'position','roi');
            loct=obj.locData.getloc({'xnm','ynm','znm','locprecnm','locprecznm'},'layer',p.targetlayer.Value,'position','roi');
            
           %determine range
           q1=0.01;
           q2=1-q1;
           minx=min(myquantilefast(locr.xnm,q1),myquantilefast(loct.xnm,q1));maxx=max(myquantilefast(locr.xnm,q2),myquantilefast(loct.xnm,q2));
           miny=min(myquantilefast(locr.ynm,q1),myquantilefast(loct.ynm,q1));maxy=max(myquantilefast(locr.ynm,q2),myquantilefast(loct.ynm,q2));
           minz=max(-5000,min(myquantilefast(locr.znm,q1),myquantilefast(loct.znm,q1)));maxz=min(5000,max(myquantilefast(locr.znm,q2),myquantilefast(loct.znm,q2)));
           
           zb=minz:p.binwidth:maxz;
           xb=minx:p.slicewidth:maxx; xb(end)=maxx;
           
           if p.setpixelsize
               pixrec=p.pixelsize;
           else
               pixrec=p.sr_pixrec;
           end
           imref=myhist2(locr.ynm,locr.xnm,pixrec,pixrec,[miny maxy],[minx maxx]);
           imtarget=myhist2(loct.ynm,loct.xnm,pixrec,pixrec,[miny maxy],[minx maxx]);
           initaxis(p.resultstabgroup,'x-y correlation');
           
           [dy,dx]=getShiftCorr(imref,imtarget,1);
           dx=dx*pixrec;
           dy=dy*pixrec;
           title([dx dy])
           ax=obj.initaxis('z-correlation');
           zpos=finddisplacementZ(locr.xnm,locr.znm,loct.xnm,loct.znm,xb,zb,5,ax);
           zpos
           
           obj.setGuiParameters(struct('foundd',[dx dy zpos]));
           clipboard('copy',sprintf([num2str(dx) ',' num2str(dy) ',' num2str(zpos) ]))

%            [~,sindr]=sort(locr.xnm);
%            [~,sindt]=sort(loct.xnm);
%            x1r=1;x1t=1;
%            ccc=zeros(1,length(zb)*2-3);
%            for k=1:length(xb)-1
%                x2r=x1r;x2t=x1t;
%                while(locr.xnm(sindr(x2r))<xb(k+1))&&x2r<length(sindr)
%                    x2r=x2r+1;
%                end
%                while(loct.xnm(sindt(x2t))<xb(k+1))&&x2t<length(sindt)
%                    x2t=x2t+1;
%                end
%         
%                zr=locr.znm(sindr(x1r:x2r-1));
%                zt=loct.znm(sindt(x1t:x2t-1));
% 
%                hr=histcounts(zr,zb);
%                ht=histcounts(zt,zb);
%                hr=hr-mean(hr);
%                 ht=ht-mean(ht);
%                ch=conv(hr,ht(end:-1:1),'full');
%                ccc=ccc+ch;
%                x1r=x2r;x1t=x2t; 
%            end
%            ax=obj.initaxis('cross-correlation');
%            zc=(-length(hr)+1:length(hr)-1)*p.binwidth;
%            plot(zc,ccc)
%            title(ax,['dz (mean): ' num2str(mean(locr.znm)-mean(loct.znm))])
%            [mc,ind]=max(ccc);
%            dh=find(ccc(ind:end)<mc/2,1,'first');
%            dh=max(3,round(dh/2));
%            zred=zc(ind-dh:ind+dh);
%            [xpos,fp]=mypeakfit(zc(ind-dh:ind+dh),ccc(ind-dh:ind+dh));
%            hold on
%            plot(zred,fp(1)*zred.^2+zred*fp(2)+fp(3),'r');
%            hold off
%            
%            title(ax,['dz (mean): ' num2str(mean(locr.znm)-mean(loct.znm)) ', cc: ' num2str(xpos)])
           
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        function makeGui(obj)
            makeGui@interfaces.DialogProcessor(obj);

        end
        function initGui(obj)
%             obj.addSynchronization('transformationfile',obj.guihandles.Tfile,'String');
        end
      
    end
end



function pard=guidef


pard.texttl.object=struct('String','Targetlayer:','Style','text');
pard.texttl.position=[2,1];
% pard.texttl.Width=0.15;

pard.targetlayer.object=struct('Style','popupmenu','String','layer1|layer2|layer3|layer4|','Value',2);
pard.targetlayer.position=[3,1];
pard.targetlayer.object.TooltipString='layer';
% pard.targetlayer.Width=0.85;

pard.texttr.object=struct('String','Referencelayer:','Style','text');
pard.texttr.position=[2,2];
% pard.texttr.Width=0.15;

pard.reflayer.object=struct('Style','popupmenu','String','layer1|layer2|layer3|layer4|');
pard.reflayer.position=[3,2];
pard.reflayer.object.TooltipString='layer';
% pard.reflayer.Width=0.85;
pard.binwidtht.object=struct('Style','text','String','Binwidth (nm)');
pard.binwidtht.position=[4,1];
pard.binwidtht.object.TooltipString='binwidth';

pard.binwidth.object=struct('Style','edit','String','10');
pard.binwidth.position=[4,2];
pard.binwidth.object.TooltipString='binwidth';

pard.slicewidtht.object=struct('Style','text','String','Slice width (nm)');
pard.slicewidtht.position=[5,1];
pard.slicewidtht.object.TooltipString='binwidth';

pard.slicewidth.object=struct('Style','edit','String','100');
pard.slicewidth.position=[5,2];
pard.slicewidth.object.TooltipString='binwidth';


pard.setpixelsize.object=struct('Style','checkbox','String','Set pixelsize (nm)','Value',0);
pard.setpixelsize.position=[1,1];
pard.setpixelsize.Width=1.5;
pard.pixelsize.object=struct('Style','edit','String','5');
pard.pixelsize.position=[1,2.5];
pard.pixelsize.Width=0.5;

pard.founddt.object=struct('Style','text','String','found dixplacemetn X,Y,Z (nm)');
pard.founddt.position=[6,1];
pard.founddt.Width=1.5;
pard.foundd.object=struct('Style','edit','String',' ');
pard.foundd.position=[6,2.5];
pard.foundd.Width=1.5;

% pard.syncParameters={{'filelist_short','dataselect',{'String'}}};
% pard.inputParameters={'currentfileinfo'};

pard.plugininfo.name='XYZ displacement';
pard.plugininfo.type='ProcessorPlugin';
end