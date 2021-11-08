classdef imageloaderTifSimple<interfaces.imageloaderSMAP
    %imageloaderMM image loader for micromanager  tiff stack files
    %   Detailed explanation goes here
    
    properties
         reader
    end
    
    methods
        function obj=imageloaderTifSimple(varargin)
            obj@interfaces.imageloaderSMAP(varargin{:});
        end
        function openi(obj,file)
        
%             obj.reader = javaObjectEDT('org.micromanager.acquisition.TaggedImageStorageMultipageTiff',fileparts(file), false, [], false, false, true);
            obj.file=file;
            warning('off','imageio:tiffmexutils:libtiffWarning');
            obj.reader=Tiff(obj.file);
%             obj.reader=bfGetReader(file);
            md=obj.getmetadata;
%             [p,f]=fileparts(file);
            obj.metadata.basefile=file;
            
            
        end
        function image=getimagei(obj,frame)
            try
                obj.reader.setDirectory(frame);
                image=obj.reader.read;
%                 image=imread(obj.file,'Index',frame);
            catch
                image=[];
            end
           
        end
        
        function closei(obj)
            obj.reader.close;
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
            %try to read out essential tiff meta data
%             allmd(end+1,:)={'Width info',obj.reader.getTag("ImageWidth")};
%             allmd(end+1,:)={'Height info',obj.reader.getTag("ImageLength")};
%             allmd(end+1,:)={'FileName',obj.file};

            obj.reader.setDirectory(1);
            imtest=obj.reader.read;
            allmd(end+1,:)={'Width info',size(imtest,2)};
            allmd(end+1,:)={'Height info',size(imtest,1)};
            allmd(end+1,:)={'FileName',obj.file};
            allmd(end+1,:)={'Roi direct',num2str([0 0 size(imtest,2) size(imtest,1)])};
%             warning('off','imageio:tiffmexutils:libtiffWarning');
%             ttt=Tiff(obj.file);
            numf=1;
            try
                desc=ttt.getTag('ImageDescription');
            catch err
                disp('no description found')
                desc='';               
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
                numf=getTifframes(obj.reader);
%                 warning('off','MATLAB:imagesci:tiffmexutils:libtiffWarning');
%                 try
%                     while ~obj.reader.lastDirectory
%                         obj.reader.nextDirectory
%                     end
%                 catch err
% %                     disp('could not determine number of images')
%                 end
%                 numf=obj.reader.currentDirectory;
            end
%             ttt.close;
            allmd(end+1,:)={'Frames',numf};
            allmd(end+1,:)={'frames direct',numf};
            obj.allmetadatatags=allmd;
                
        
        end
        
    end
    
end

function numf=getTifframes(reader)
try
    reader.setDirectory(2^16);
catch err
end
numf=reader.currentDirectory-1;

% if nargin<2
%     df=1;
% end
% currentdir=reader.currentDirectory;
% try
%     while ~reader.lastDirectory
% %         if df>1
%             reader.setDirectory(reader.currentDirectory+df);
%             newdir=reader.currentDirectory;
%             if newdir<currentdir
%                 break
%             end
%             currentdir=newdir;
% %         else
% %             reader.nextDirectory
% %         end
%     end
% catch err
%     err
% %                     disp('could not determine number of images')
% end
% numf=reader.currentDirectory;
% if df>1
%     numf=getTifframes(reader,ceil(df/10));
% end

end
