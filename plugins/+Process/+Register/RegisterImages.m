classdef RegisterImages<interfaces.DialogProcessor
    properties
        transformation
    end
    methods
        function obj=RegisterImages(varargin)        
            obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'sr_pixrec','sr_pos','sr_size'};
            obj.showresults=true;
            obj.undo=true;
        end
        
        function out=run(obj,p)
            out=[];
            targetlayer=p.targetselect.Value;
            reflayer=p.refselect.Value;
            ll=length(obj.locData.layer);
            if targetlayer>ll||reflayer>ll
                out.error=('selected layer does not exist or has not been calculated');
                return;
            end
            
            if p.setpixelsize
                locsref=obj.locData.getloc({'xnm','ynm','locprecnm'},'layer',reflayer);
                locstarget=obj.locData.getloc({'xnm','ynm','locprecnm'},'layer',targetlayer);
                xr=[p.sr_pos(1)-p.sr_size(1) p.sr_pos(1)+p.sr_size(1)];
                yr=[p.sr_pos(2)-p.sr_size(2) p.sr_pos(2)+p.sr_size(2)];
                imref=myhist2(locsref.ynm,locsref.xnm,p.pixelsize,p.pixelsize,yr,xr);
                imtarget=myhist2(locstarget.ynm,locstarget.xnm,p.pixelsize,p.pixelsize,yr,xr);
                pixrec=p.pixelsize;
            else
                pixrec=p.sr_pixrec;
                imref=sum(obj.locData.layer(reflayer).images.finalImages.image,3);
                imtarget=sum(obj.locData.layer(targetlayer).images.finalImages.image,3);
            end
            if sum(imref(:))==0||sum(imtarget(:))==0
                out.error=('reconstructed images are empty.');
                return;
            end
            initaxis(p.resultstabgroup,'correlation');
             [dy,dx]=getShiftCorr(imref,imtarget,1);
             dx=dx*pixrec;
             dy=dy*pixrec;
                dxt=obj.getPar('shiftxy_min','layer',targetlayer);
                dyt=obj.getPar('shiftxy_max','layer',targetlayer);
%             dxr=obj.getPar('shiftxy_min','layer',reflayer);
%             dyr=obj.getPar('shiftxy_max','layer',reflayer);
            if p.correctcoordinates
                
                tch=obj.locData.loc.channel==p.targetchannels(1);
                for k=2:length(p.targetchannels)
                    tch=tch|(obj.locData.loc.channel==p.targetchannels(k));
                end
                tfile=p.dataselect.Value==obj.locData.loc.filenumber;
                tch=tch&tfile;
                obj.locData.loc.xnm(tch)=obj.locData.loc.xnm(tch)+dx;
                obj.locData.loc.ynm(tch)=obj.locData.loc.ynm(tch)+dy;
                obj.locData.regroup;
            else

                obj.setPar('shiftxy_min','layer',targetlayer,dx+dxt(1))
                obj.setPar('shiftxy_max','layer',targetlayer,dy+dyt)
            end
            notify(obj.P,'sr_render');
            
            %make locTransform
            A=[1,0,0;0,1,0;(dx+dxt(1))/1000,(dy+dyt)/1000,1];
            trafo=interfaces.LocTransform;
            trafo.makeAffine2d(A);
            obj.transformation=trafo;
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function initGui(obj)
            obj.addSynchronization('transformationfile',obj.guihandles.Tfile,'String');
        end
        function savebutton(obj,a,b)
            fn=obj.guihandles.Tfile.String;
            [f,path]=uiputfile(fn,'Select transformation file _T.mat');
            if f
                obj.guihandles.Tfile.String=[path f];
                transformation=obj.transformation;
                save([path,f],'transformation');
            end      
        end 
    end
end




function pard=guidef(obj)
pard.texta.object=struct('String','target','Style','text');
pard.texta.position=[1,1];
pard.textb.object=struct('String','reference','Style','text');
pard.textb.position=[1,2];
pard.targetselect.object=struct('Style','popupmenu','String','layer1|layer2|layer3|layer4|layer5','Value',2);
pard.targetselect.position=[2,1];
pard.targetselect.load=false;
pard.refselect.object=struct('Style','popupmenu','String','layer1|layer2|layer3|layer4|layer5','Value',1);
pard.refselect.position=[2,2];
pard.refselect.load=false;

pard.correctcoordinates.object=struct('Style','checkbox','String','Correct coordinates','Value',0);
pard.correctcoordinates.position=[3,1];
pard.correctcoordinates.Width=2;

pard.targetchannelst.object=struct('Style','text','String','Channels:');
pard.targetchannelst.position=[4,1];
pard.targetchannelst.Width=.5;

pard.targetchannels.object=struct('Style','edit','String','1 ');
pard.targetchannels.position=[4,1.5];
pard.targetchannels.Width=.5;

pard.setpixelsize.object=struct('Style','checkbox','String','Set pixelsize (nm)','Value',0);
pard.setpixelsize.position=[1,3];
pard.setpixelsize.Width=1.5;
pard.pixelsize.object=struct('Style','edit','String','5');
pard.pixelsize.position=[1,4.5];
pard.pixelsize.Width=0.5;

pard.dataselect.object=struct('Style','popupmenu','String','File');
pard.dataselect.position=[5,1];
pard.dataselect.object.TooltipString='choose localization file data set';

pard.syncParameters={{'filelist_short','dataselect',{'String'}}};


pard.Tfile.object=struct('Style','edit','String','settings/temp/temp_T.mat');
pard.Tfile.position=[8,1];
pard.Tfile.Width=3;

pard.savebutton.object=struct('Style','pushbutton','String','save T','Callback',@obj.savebutton);
pard.savebutton.position=[8,4];

pard.plugininfo.description=sprintf(['Register Images calculates shift between rendered images in two layers and writes this shift into the channel tab of the target layer.'...
    '\n Adjust pixel size and FoV so that the shift can be calculated from the reconstructed image.'...
    '\n Transformation can also be saved for later use with Apply Transform']);
pard.plugininfo.name='Register Images';
pard.plugininfo.type='ProcessorPlugin';
end