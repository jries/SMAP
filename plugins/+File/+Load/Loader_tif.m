classdef Loader_tif<interfaces.DialogProcessor
%     Loads diffrction-limited tiff images and associates them to an
%     already loaded SMAP data set.
    methods
        function obj=Loader_tif(varargin)        
                obj@interfaces.DialogProcessor(varargin{:}) ;
                obj.inputParameters={'mainGui','filelist_short'};
        end
        
        function out=load(obj,p,file,mode)
            if nargin<4
                mode=getfilemode(file);
            end
            p=obj.getAllParameters;
            loadfile(obj,p,file,mode);
        end
        
        function pard=guidef(obj)
            pard=guidef;
        end
        function run(obj,p)
            [f,path]=uigetfile(obj.info.extensions);
            if exist([path f],'file')
                obj.load(p,[path f]);
                initGuiAfterLoad(obj);
                out.file=[f,path];
            else
                out.error='file not found. Cannot be loaded.';
            end
        end
        function clear(obj,file,isadd)
                obj.locData.clear('filter');
        end        

    end
end




function pard=guidef
info.name='tif loader';
info.extensions={'*.tif';'*.*'};
info.dialogtitle='select any Tif file';

pard.stack.object=struct('Style','checkbox','String','load stack');
pard.stack.position=[1,1];
pard.stack.Width=1.3;
pard.stack.TooltipString='Load image stack (otherwise load single image)';

pard.mirrortif.object=struct('Style','checkbox','String','Mirror (for EM gain)');
pard.mirrortif.position=[1,2.3];
pard.mirrortif.Width=1.3;
pard.mirrortif.TooltipString='Mirror image. Might be required if EM gain was used';

pard.autocontrast.object=struct('Style','checkbox','String','auto contrast');
pard.autocontrast.position=[1,3.6];
pard.autocontrast.Width=1.3;

pard.abberior.object=struct('Style','checkbox','String','Abberior. Offset x,y: ');
pard.abberior.position=[2,1];
pard.abberior.Width=2;

pard.offsetsx.object=struct('Style','edit','String','0');
pard.offsetsx.position=[2,3];
pard.offsetsx.Width=1;
pard.offsetsy.object=struct('Style','edit','String','0');
pard.offsetsy.position=[2,4];
pard.offsetsy.Width=1;


pard.plugininfo=info;
pard.plugininfo.type='LoaderPlugin';
pard.plugininfo.description='Loads diffrction-limited tiff images and associates them to an already loaded SMAP data set.';
end

function loadfile(obj,p,file,mode)  
if ~strcmp(mode,'tif')
   disp('no recognized Tiff file')
   return
end

if obj.locData.files.filenumberEnd==0
    [path,f]=fileparts(file);
    filename=[path filesep f];
   obj.locData.addfile(filename) 
end
if length(obj.locData.files.file)>1
    s=p.filelist_short.String;
    f=listdlg('ListString',s);
else
    f=1;
end

if isfield(p,'stack')&&p.stack
    images=chooseTifImage(file,obj.P);
else
    images=gettif(file);
end

if isfield(p,'mirrortif' ) && p.mirrortif
    for k=1:length(images)
        images(k).image=images(k).image(:,end:-1:1);
        roih=images(k).info.roi;
        % images(k).info.roi(1)=512-roih(1)-roih(3);
        % added for Sarah for M2. MIght not be the same for other cameras.
    end
end
numimages=length(images);
if isfield(obj.locData.files.file(f),'numberOfTif')
    tiffold=obj.locData.files.file(f).numberOfTif;
else
    tiffold=0;
end
obj.locData.files.file(f).numberOfTif=tiffold+numimages;
% imout=gettif(file);
if p.abberior
    imf=imfinfo(file);
    px=1/imf.XResolution;
    images(1).info.roi=[p.offsetsx/px, p.offsetsy/px, size(images.image,1), size(images.image,2)];
    images(1).info.cam_pixelsize_um=px;
else
    for k=1:numimages
        images(k).info.cam_pixelsize_um=obj.locData.files.file(f).info.cam_pixelsize_um;
    end
end
if p.autocontrast
    image=double(images.image);
    image=image-min(image);
    image=image/max(image(:));
    images.image=image;
end

obj.locData.files.file(f).tif(tiffold+1:obj.locData.files.file(f).numberOfTif)=images;

  initGuiAfterLoad(obj)    
end
        


