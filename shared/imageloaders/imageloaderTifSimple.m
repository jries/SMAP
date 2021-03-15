classdef imageloaderTifSimple<interfaces.imageloaderSMAP
    %imageloaderMM image loader for micromanager  tiff stack files
    %   Detailed explanation goes here
    
    properties
%         reader
    end
    
    methods
        function obj=imageloaderTifSimple(varargin)
            obj@interfaces.imageloaderSMAP(varargin{:});
        end
        function openi(obj,file)
        
%             obj.reader = javaObjectEDT('org.micromanager.acquisition.TaggedImageStorageMultipageTiff',fileparts(file), false, [], false, false, true);
            obj.file=file;
%             obj.reader=bfGetReader(file);
            md=obj.getmetadata;
%             [p,f]=fileparts(file);
            obj.metadata.basefile=file;
            
        end
        function image=getimagei(obj,frame)
            try
            image=imread(obj.file,'Index',frame);
            catch
                image=[];
            end
           
        end
        
        function closei(obj)
%             obj.reader.close
%             clear(obj.reader)
        end
        
        function image=getimageonline(obj,number)
            image=obj.getimage(number);
            if isempty(image)&&obj.onlineAnalysis 
                    disp('wait')
%                     obj.reader.close;
%                     delete(obj.reader)
                    pause(obj.waittime*2)
%                     obj.reader = javaObjectEDT('org.micromanager.acquisition.TaggedImageStorageMultipageTiff',fileparts(obj.file), false, [], false, false, true);
                    image=obj.getimage(number);
            end
        end
        
        function allmd=getmetadatatagsi(obj)
            allmd={'Format','SimpleTif'};
            imtest=imread(obj.file,'Index',1);
            allmd(end+1,:)={'Width info',size(imtest,2)};
            allmd(end+1,:)={'Height info',size(imtest,1)};
            allmd(end+1,:)={'FileName',obj.file};
            warning('off','imageio:tiffmexutils:libtiffWarning');
            ttt=Tiff(obj.file);
            numf=1;
            try
            desc=ttt.getTag('ImageDescription');
            catch err
                disp('no description found')
                desc='';
                %manually find frames
                
%                 ttt.setDirectory(numf)
%                 while ~ttt.lastDirectory
%                     numf=numf+1;
%                     ttt.setDirectory(numf);
%                 end
                    
            end
           
            ind=strfind(desc,'images=');
            
            if ~isempty(ind)
                numf=sscanf(desc(ind:end),'images=%f');
            else
                ind=strfind(desc,'slices=');
                if ~isempty(ind)
                    numf=sscanf(desc(ind:end),'slices=%f');
                end
            end
            if numf==1 %might not have worked
                disp('could not determine number of images. Directly evaluate file, might take time.');
                warning('off','MATLAB:imagesci:tiffmexutils:libtiffWarning');
                try
                    while ~ttt.lastDirectory
                        ttt.nextDirectory
                    end
                catch err
                    disp('could not determine number of images')
                end
                numf=ttt.currentDirectory;
            end
            ttt.close;
            allmd(end+1,:)={'Frames',numf};
            allmd(end+1,:)={'frames direct',numf};
            obj.allmetadatatags=allmd;
                
        
        end
        
    end
    
end


