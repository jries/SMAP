classdef Loader_tif<interfaces.DialogProcessor
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
            [f,p]=uigetfile(obj.info.extensions);
            obj.load(p,[p f]);
            initGuiAfterLoad(obj);
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
pard.stack.Width=2;
pard.stack.TooltipString='Load image stack (otherwise load single image)';

pard.mirrortif.object=struct('Style','checkbox','String','Mirror (for EM gain)');
pard.mirrortif.position=[1,3];
pard.mirrortif.Width=2;
pard.mirrortif.TooltipString='Mirror image. Might be required if EM gain was used';

pard.plugininfo=info;
pard.plugininfo.type='LoaderPlugin';
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
        images(k).info.roi(1)=512-roih(1)-roih(3);
    end
end
numimages=length(images);
tiffold=obj.locData.files.file(f).numberOfTif;
obj.locData.files.file(f).numberOfTif=tiffold+numimages;
% imout=gettif(file);
for k=1:numimages
    images(k).info.cam_pixelsize_um=obj.locData.files.file(f).info.cam_pixelsize_um;
end

obj.locData.files.file(f).tif(tiffold+1:obj.locData.files.file(f).numberOfTif)=images;

  initGuiAfterLoad(obj)    
end
        


