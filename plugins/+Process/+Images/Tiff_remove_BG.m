classdef Tiff_remove_BG<interfaces.DialogProcessor
    methods
        function obj=Tiff_remove_BG(varargin)      
            obj@interfaces.DialogProcessor(varargin{:}) ;  
        end
        
        function out=run(obj,p)
            out=[];
            obj.setPar('undoModule','Tiff_remove_BG');
            notify(obj.P,'backup4undo');
            file=p.dataselect.Value;
            tiffn=p.tiffselect.Value;
            tiff=obj.locData.files.file(file).tif(tiffn);
            img=tiff.image;
            bg=mywaveletfilter(double(img),p.waveletlevel,true);
            imgc=double(img)-bg;
            [path,fileh, ext]=fileparts(tiff.info.name);
            fileh=['BGc_' fileh ext];
            tiff.info.name=fullfile(path,fileh);
            tiff.image=imgc;
            
            if ~p.previewcheck
                if p.overwrite
                    obj.locData.files.file(file).tif(tiffn)=tiff;
                else
                    obj.locData.files.file(file).tif(end+1)=tiff;
                end
            end
            ax=obj.initaxis('image');
            imagesc(img,'Parent',ax);
            ax3=obj.initaxis('bg');
            imagesc(bg,'Parent',ax3);
            ax2=obj.initaxis('image-bg');
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

pard.waveletlevelt.object=struct('Style','text','String','wavelet level: ');
pard.waveletlevelt.position=[3,1];
pard.waveletlevel.object=struct('Style','edit','String','3');
pard.waveletlevel.position=[3,2];
pard.waveletlevel.object.TooltipString='Large value: large scale filtering. Typical 2-4';
pard.previewcheck.object=struct('Style','checkbox','String','preview','Value',1);
pard.previewcheck.position=[4,1];
pard.previewcheck.Width=2;
pard.previewcheck.object.TooltipString='Preview if checked. Otherwise the image is saved.';
pard.overwrite.object=struct('Style','checkbox','String','replace, not append');
pard.overwrite.position=[5,1];
pard.overwrite.Width=2;
pard.overwrite.object.TooltipString='If checked, the original image is replaced. Otherwise the BG corrected image is appended to the tiffs.';

pard.syncParameters={{'filelist_short','dataselect',{'String'}}};
pard.plugininfo.type='ProcessorPlugin';
end