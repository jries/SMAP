classdef Rescale_tiff<interfaces.DialogProcessor
%     Change position and magnification of an associated Tiffi mage
    methods
        function obj=Rescale_tiff(varargin)      
            obj@interfaces.DialogProcessor(varargin{:}) ;  
        end
        
        function out=run(obj,p)
            out=[];
            obj.setPar('undoModule','Rescale_tiff');
            notify(obj.P,'backup4undo');
            file=p.dataselect.Value;
            tiffn=p.tiffselect.Value;
            tiff=obj.locData.files.file(file).tif(tiffn);
            img=tiff.image;
            Mx=p.magnification(1);
            My=p.magnification(end);
            pixelsize=tiff.info.cam_pixelsize_um*1000;
            sx=p.shift(1)/1000;
            sy=p.shift(end)/1000;
            a=p.rot;
            
            M=[Mx*cosd(a) -My*sind(a) 0
              Mx*sind(a) My*cosd(a) 0
                sx sy 1];
            tform=affine2d(M);
            imgc=transformImage(tform,img,pixelsize,tiff.info.roi);
            
            [path,fileh, ext]=fileparts(tiff.info.name);
            fileh=['T_' fileh ext];
            tiff.info.name=fullfile(path,fileh);
            tiff.image=imgc;
            
            if 1 %~p.previewcheck
                if p.overwrite
                    obj.locData.files.file(file).tif(tiffn)=tiff;
                else
                    obj.locData.files.file(file).tif(end+1)=tiff;
                end
            end
            ax=obj.initaxis('image');
            imagesc(img,'Parent',ax);
            ax2=obj.initaxis('image transformed');
            imagesc(imgc,'Parent',ax2);            
             
            
             
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function dataselect_callback(obj,object,b)
            tifs=obj.locData.files.file(object.Value).tif;
            info=[tifs.info];
            for k=1:length(info)
                [~,tifnames{k}]=fileparts(info(k).name);
            end
           
            obj.guihandles.tiffselect.String=tifnames;
            obj.guihandles.tiffselect.Value=max(1,min(obj.guihandles.tiffselect.Value,length(tifnames)));
            
        end
    end
end

function pard=guidef(obj)
pard.dataselect.object=struct('Style','popupmenu','String','File','Callback',{{@obj.dataselect_callback}});
pard.dataselect.position=[1,1];
pard.dataselect.object.TooltipString='file for target localizations';
pard.dataselect.Width=2;
pard.tiffselect.object=struct('Style','popupmenu','String','empty');
pard.tiffselect.position=[2,1];
pard.tiffselect.object.TooltipString='select Tiff. Use File selector first to populate this menu.';
pard.tiffselect.Width=2;

pard.magnificationt.object=struct('Style','text','String','Magnification Mx My');
pard.magnificationt.position=[3,1];
pard.magnification.object=struct('Style','edit','String','1 1');
pard.magnification.position=[3,2];
pard.magnification.object.TooltipString='Magnification in x and y';

pard.shiftt.object=struct('Style','text','String','shift x y (nm)');
pard.shiftt.position=[4,1];
pard.shift.object=struct('Style','edit','String','0 0');
pard.shift.position=[4,2];
pard.shift.object.TooltipString='Shift in x and y (nm)';

pard.rott.object=struct('Style','text','String','rotation angle (deg)');
pard.rott.position=[5,1];
pard.rot.object=struct('Style','edit','String','0');
pard.rot.position=[5,2];
pard.rot.object.TooltipString='Rotation angle in deg';


% pard.previewcheck.object=struct('Style','checkbox','String','preview','Value',1);
% pard.previewcheck.position=[6,1];
% pard.previewcheck.Width=2;
% pard.previewcheck.object.TooltipString='Preview if checked. Otherwise the image is saved.';
pard.overwrite.object=struct('Style','checkbox','String','replace, not append');
pard.overwrite.position=[7,1];
pard.overwrite.Width=2;
pard.overwrite.object.TooltipString='If checked, the original image is replaced. Otherwise the BG corrected image is appended to the tiffs.';

% pard.inputParameters={'cam_pixnm'
pard.syncParameters={{'filelist_short','dataselect',{'String'}}};
pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.description='Change position and magnification of an associated Tiffi mage';
end