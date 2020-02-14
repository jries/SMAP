classdef SetTiffLocation<interfaces.DialogProcessor
% Sets location of associated raw data tiff files. Useful for batch
% processing.
    properties 
    end
    methods
        function obj=SetTiffLocation(varargin)   
            obj@interfaces.DialogProcessor(varargin{:}) ;
        end
        
        function out=run(obj,p)
            out=[];
            tiffile=getrawtifpath(obj.locData);
            [file,pfad]=uigetfile(tiffile);
            if file
                obj.locData.files.file(1).info.imagefile=[pfad file];
                if p.savesmlfile
                    obj.locData.savelocs(obj.locData.files.file(1).name);
                end
            end
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end




function pard=guidef(obj)
pard.savesmlfile.object=struct('Style','checkbox','String','save SML file','Value',1);
pard.savesmlfile.position=[1,1];
pard.savesmlfile.Width=4;



pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.description='Sets location of associated raw data tiff files. Useful for batch processing.';
end



function tiffile=getrawtifpath(locData)
    if isfield(locData.files.file(1).info,'imagefile')
        tiffile=locData.files.file(1).info.imagefile;
    else
        tiffile=locData.files.file(1).info.basefile;
    end
    if ~exist(tiffile,'file')
        disp('Get2CIntImagesWF ine 40: check if it works')
        tiffile=strrep(tiffile,'\','/');
        ind=strfind(tiffile,'/');
        for k=1:length(ind)
            tiffileh=[tiffile(1:ind(k)) '_b_' tiffile(ind(k)+1:end)];
            if exist(tiffileh,'file')
                tiffile=tiffileh;
            end
        end
    end
    if ~exist(tiffile,'file')
        tiffile=strrep(locData.files.file(1).name,'_sml.mat','.tif');
    end
end
