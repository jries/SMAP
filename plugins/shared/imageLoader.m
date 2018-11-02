classdef imageLoader<handle
    properties
        file
        imageMode
        onlineAnalysis=0;
        waittime=1;
        separate
        separatefiles
        currentImageNumber
        info
        stack
    end
    methods
        function obj=imageLoader(varargin)
            if nargin>0
                obj.setFile(varargin{1});
            end
        end
        function setFile(obj,file)
            
            obj.info=getimageinfo(file);
            obj.file=obj.info.filename;
            
            switch obj.info.format
                case 'stackTif'
                    obj.imageMode=1;
                    obj.stack.currentfile=1;
                    obj.updatestack(obj.info);
%                     sd=warning;
                    warning('off','MATLAB:imagesci:tiffmexutils:libtiffWarning')
                    obj.stack.tiffh=Tiff(obj.file,'r');
%                     warning(sd);
                    obj.currentImageNumber=0;
%                     warning('off','MATLAB:imagesci:tiffmexutils:libtiffWarning');
                case 'separateTif'
                    obj.imageMode=0;
                    obj.separatefiles=obj.info.files;
                    obj.separate.numfiles=length(obj.info.files);
                    obj.separate.path=obj.info.path;
                    obj.currentImageNumber=0;
                    obj.separate.fmt_s=imformats('tif');
            end
%             obj.ftiff=Tiff(obj.file);
        end
        function updatestack(obj,info)
%             obj.stack.tiffh=info.tiffh;
%             obj.stack.lastread=false;
            numframes=[info.allfiles(1:end-1).numberOfFrames];
            obj.stack.imageoffset=cumsum([0 numframes]);
            obj.stack.lastimage=cumsum([info.allfiles(1:end).numberOfFrames]);
            obj.stack.files={info.allfiles(:).name};
        end
        function info=getinfo(obj)
            info=obj.info;           
        end
        function image=readNext(obj)
            obj.currentImageNumber=obj.currentImageNumber+1;
            if obj.imageMode==0 %separate Tiff
                image=obj.readSeparate(obj.currentImageNumber);  
            elseif obj.imageMode==1 %tiffstack
                image=obj.readstack(obj.currentImageNumber);     
            end
        end
        function image=readSeparate(obj,number)
                image=readseparate(obj,number);           
        end
        function [imagenumber,filenumberc,filenumber]=getstacknumber(obj,frame)
            filenumber=find(frame<=obj.stack.lastimage,1,'first');
            filenumberc=filenumber;
            if isempty(filenumber)
                filenumberc=length(obj.stack.lastimage);
            end
            imagenumber=frame-obj.stack.imageoffset(filenumberc);
        end
        
        function image=readstack(obj,imagenumber,recursions)
            if nargin<3
                recursions=0;
            end
            maxrecursions=1;
            image=[]; %default, if nothing can be read: end.

            [imagenumber,filenumber]=obj.getstacknumber(imagenumber);
            if filenumber~=obj.stack.currentfile
                obj.stack.tiffh.close;
                newfile=obj.stack.files{filenumber};
                obj.info=getimageinfo(newfile);
                obj.updatestack(obj.info);
                obj.stack.tiffh=Tiff(newfile,'r');
                obj.stack.currentfile=filenumber;
            end
            try
                th=obj.stack.tiffh;
                th.setDirectory(imagenumber);
                image=th.read;                    
            catch
                
                if obj.onlineAnalysis
                    if recursions<maxrecursions
                        pause(obj.waittime*2);
                        obj.stack.currentfile=-1; %update stuff
                        image=obj.readstack(imagenumber,recursions+1);
                    else
                        disp('after waiting no more files')
                    end
                else
                    disp('image out of range');
                end
            end
        end
        
        function image=getImage(obj,number)
            obj.setImageNumber(number-1);
            image=obj.readNext;
        end
        
        function setImageNumber(obj,number)
            obj.currentImageNumber=number;
%              if obj.imageMode==1
%                  try
%                      obj.stack.tiffh.setDirectory(number+1)
%                      obj.stack.lastread=false;
%                  catch
%                      disp('number outside tiff range, set to 1')
%                      obj.setImageNumber(1);
%                  end
%              end
        end
    end
end


function image=readseparate(obj,number)
separate=obj.separate;
% lenfiles=length(separate.files);
lenfiles=separate.numfiles;
if lenfiles<number 
    if obj.onlineAnalysis %ask for image that is not in list
        lastfile= obj.separatefiles{lenfiles};
        thisname= generateFileName(lastfile,lenfiles,obj.info.numberNameRange,number);
        thisfile=[separate.path filesep thisname];
        if ~exist(thisfile,'file')
            disp('wait')
            pause(obj.waittime*2)
        end
        if ~exist(thisfile,'file')
            image=[];
            return
        else
            if number>lenfiles
                obj.separatefiles{number+1000}='';
            end
            obj.separatefiles{number}=thisname;
            obj.separate.numfiles=max(lenfiles,number);
            obj.info.numberOfFrames=max(lenfiles,number);
        end 
    else
        image=[];
        return
    end
else
    thisfile=[separate.path filesep obj.separatefiles{number}];
end
image=myimread(thisfile,separate.fmt_s);     
end


function newfile=generateFileName(oldfile, oldfilenumber,indfbar,newfilenumber)
oldfilenamenumber=str2double(oldfile(indfbar(1):indfbar(2)));
newfilenamenumber=oldfilenamenumber-oldfilenumber+newfilenumber;
lenfield=indfbar(2)-indfbar(1)+1;
newfilestr=num2str(newfilenamenumber,['%0' num2str(lenfield) 'i']);
newfile=[oldfile(1:indfbar(1)-1) newfilestr oldfile(indfbar(2)+1:end)];
end