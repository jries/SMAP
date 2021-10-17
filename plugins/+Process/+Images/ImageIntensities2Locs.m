classdef ImageIntensities2Locs<interfaces.DialogProcessor
%     Change position and magnification of an associated Tiffi mage
    methods
        function obj=ImageIntensities2Locs(varargin)      
            obj@interfaces.DialogProcessor(varargin{:}) ;  
        end
        
        function out=run(obj,p)
            out=[];
            obj.setPar('undoModule','Rescale_tiff');
            notify(obj.P,'backup4undo');
            locs=obj.locData.getloc({'xnm','ynm','filenumber','tiffintensity'},'grouping','ungrouped');
            %calculate in tiff pixels
            file=p.dataselect.Value;
            tiff=p.tiffselect.Value;
            infile=locs.filenumber==file;
            tiffinfo=obj.locData.files.file(file).tif(tiff).info;
             img=obj.locData.files.file(file).tif(tiff).image;
            xpix=locs.xnm(infile)/tiffinfo.cam_pixelsize_um(1)/1000-tiffinfo.roi(1); %check if reall +1
            ypix=locs.ynm(infile)/tiffinfo.cam_pixelsize_um(2)/1000-tiffinfo.roi(2);
            badx=xpix<1 | xpix>size(img,2); %or 2?
            bady=ypix<1 | ypix>size(img,1); %or 1?
            xpix(badx)=1;
            ypix(bady)=1;
            linind=sub2ind(size(img),round(ypix),round(xpix));
            
            intensity=img(linind);
            intensity(badx | bady)=0;
            if isempty(locs.tiffintensity)
                outloc=zeros(size(locs.xnm),'single');
            else
                outloc=locs.tiffintensity;
            end
            outloc(infile)=intensity;
            obj.locData.setloc('tiffintensity',outloc)
            obj.locData.regroup;
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



% pard.inputParameters={'cam_pixnm'
pard.syncParameters={{'filelist_short','dataselect',{'String'}}};
pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.description='Write intensity of tiff image to localization data';
end