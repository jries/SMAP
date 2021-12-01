classdef imageloaderSMAP<interfaces.GuiParameterInterface
    %imageloaderSMAP Superclass for image loader classes
    
    properties
        metadata=interfaces.metadataSMAP;
        file='';
        onlineAnalysis=false;
        waittime=5;
        currentImageNumber;
        allmetadatatags;
        calibrationFile;
        allowmultiplefiles=true;
        ismultichannel=false;
        multiloader={};
        multiloadermetadata;
    end
    methods
       function obj=imageloaderSMAP(varargin)
           obj.metadata=interfaces.metadataSMAP;
       
           if nargin>2 && ~isempty(varargin{3})
               if isa(varargin{3},'interfaces.ParameterData')
                    obj.P=varargin{3};
                    obj.calibrationFile=obj.getGlobalSetting('cameraSettingsFile');
               else
                   obj.calibrationFile=varargin{3};
               end
           else
                   obj.calibrationFile=[obj.getPar('SettingsDirectory') '/cameras.mat']; 
           end
           if nargin>1 && ~isempty(varargin{2})
                obj.updatemetadata(varargin{2});
                obj.multiloadermetadata=varargin{2};
           end
      
            if nargin>0 && ~isempty(varargin{1})
                obj.open(varargin{1});
            end

        end
    end
    methods (Abstract=true)
        openi(obj,filename);
        image=getimagei(obj,frame);
        allmd=getmetadatatagsi(obj)
    end
    methods
        function image=getimage(obj,frame)
            if obj.ismultichannel
                for k=length(obj.multiloader):-1:1
                    image{k}=obj.multiloader{k}.getimagei(frame);
                end
                image=combineimages(image);
            else
                image=obj.getimagei(frame);
            end
        end
        function image=readNext(obj)
           obj.currentImageNumber=obj.currentImageNumber+1;
           if obj.ismultichannel
                for k=length(obj.multiloader):-1:1
                    image{k}=obj.multiloader{k}.getimageonline(obj.currentImageNumber);
                end
                image=combineimages(image);
           else
                image=obj.getimageonline(obj.currentImageNumber);
           end
        end
        
        function image=getimageonline(obj,number)
            
           if obj.ismultichannel
                for k=length(obj.multiloader):-1:1
                    image{k}=obj.getimage(number);
                    if isempty(image{k})&&obj.onlineAnalysis 
                            disp('wait')
                            pause(obj.waittime*2)
                            image{k}=obj.getimage(number);
                    end                       
                end
                image=combineimages(image);
           else
                image=obj.getimage(number);
                if isempty(image)&&obj.onlineAnalysis 
                        disp('wait')
                        pause(obj.waittime*2)
                        image=obj.getimage(number);
                end
           end
        end
        
        function setImageNumber(obj,number)
            obj.currentImageNumber=number;
        end
        function images=getmanyimages(obj,varargin)
            if obj.ismultichannel
                for k=length(obj.multiloader):-1:1
                    images{k}=obj.multiloader{k}.getmanyimagesi(varargin{:});
                end
                images=combineimages(images);
            else
                images=obj.getmanyimagesi(varargin{:});
            end
                
        end
        function images=getmanyimagesi(obj,numbers,format)
            loadfun=@obj.getimageonline;
            if nargin<3
                format='cell';
            end
            if nargin<2|| isempty(numbers)
                numbers=1:obj.metadata.numberOfFrames;
            end
            switch format
                case 'cell'
                    for k=length(numbers):-1:1
                        images{k}=loadfun(numbers(k));
                        if isempty(images{k})
                            loadfun=@obj.getimage;
                        end
                    end
                case 'mat'
                    for k=length(numbers):-1:1
                        imh=loadfun(numbers(k));
                        if ~isempty(imh)
                        images(:,:,k)=imh;
                        else
                            loadfun=@obj.getimage;
                        end
                    end
            end
            
        end
        function file=getbasefile(obj)
            file=obj.file;
        end
        function updatemetadata(obj, md)
            if isempty(md)
                return;
            end
            fn=fieldnames(md);fn2=properties(obj.metadata);
           fna=intersect(fn,fn2);
           obj.metadata=copyfields(obj.metadata,md,fna);
        end
        
        function val=gettag(obj,tag)
            if isempty(obj.allmetadatatags)
                obj.allmetadatatags=obj.getmetadatatags;
            end
            val=[];
            if ~isempty(obj.allmetadatatags)
                ind=find(strcmp(obj.allmetadatatags(:,1),tag),1,'first');
                if ~isempty(ind)
                    val=obj.allmetadatatags{ind,2};
                end
            end
        end
        
        function metao=getmetadata(obj)
            getmetadatacam(obj);
             obj.metadata.basefile=obj.getbasefile;
             metao=obj.metadata;   
        end
        function metao=getmetadatacam(obj)
            try
                camfile=obj.getGlobalSetting('cameraSettingsFile');
            catch err
                camfile=obj.calibrationFile;
                display('could not find camera file in global settings. Using default file.')
            end
%             try
                usedef=obj.getPar('useDefaultCam');
%             catch
%                 usedef=true;
%             end
            if isempty(usedef)
                usedef=true;
            end
            [md,~,~,error]=getCameraCalibration(obj,[],usedef,findsettingsfile(camfile,obj));
            if ~isempty(error)
                obj.setPar('errorindicator',error)
            else
                obj.setPar('errorindicator','clear') %when it works, clear message
            end
            if isempty(md)
                metao=[];
                return
            end
            fn=fieldnames(md);
            for k=1:length(fn)
                if ~isempty(md.(fn{k}))&&isprop(obj.metadata,fn{k})
                    obj.metadata.(fn{k})=md.(fn{k});
                    obj.metadata.assigned.(fn{k})=true;
                end
            end
            obj.metadata.allmetadata=obj.allmetadatatags;
            obj.metadata.imagefile=obj.file;
%             if isempty(obj.metadata.basefile)
               
%             end
        metao=obj.metadata;
        end
        function open(obj,file)
           if contains(file,';')
               file=strsplit(file{1},';');
           end
            if iscell(file) %multiple channels in multiple files
                obj.ismultichannel=true;
                for k=1:length(file)
                    obj.multiloader{k}= obj.newimageloader(file{k});
                end
                obj.metadata=obj.multiloader{1}.metadata; %copy metadata and infor from first
                obj.file=file;
            else
                %if '_q'
                [path,f,ext]=fileparts(file);
                indq=strfind(f,'_q');
                if ~isempty(indq)&&obj.allowmultiplefiles
                    allfiles=dir([path filesep f(1:indq+1)  '*' ext]);
                    for k=1:length(allfiles)
                        files{k}=[path filesep allfiles(k).name];
%                         disp(allfiles(k).name(indq:end))
                    end
                    files=sort(files);
                    obj.open(files);
                else
                    obj.openi(file)
                end
            end
            
        end
        function close(obj)
            if obj.ismultichannel
                for k=1:length(obj.multiloader)
                    obj.multiloader{k}.closei;
                end
            else
                obj.closei;
            end
        end
        function closei(obj)
            display(['close not implemented in ' class(obj)])
        end
        function allmd=getmetadatatags(obj)
            if obj.ismultichannel
                for k=1:length(obj.multiloader)
                    allmd=obj.multiloader{k}.getmetadatatags;
                end
                obj.metadata=obj.multiloader{1}.metadata;
            else
                allmd=obj.getmetadatatagsi;
            end
        end
        
        function il=newimageloader(obj,file)
            il=eval(class(obj));
            il.P=obj.P;
            il.calibrationFile=obj.calibrationFile;
            il.onlineAnalysis=obj.onlineAnalysis;
            il.waittime=obj.waittime;
            if ~isempty(obj.multiloadermetadata) 
                il.updatemetadata(obj.multiloadermetadata)
            end
            il.allowmultiplefiles=false;
            il.open(file);
        end
        function prefit(obj)
        end
    end
    
end


function imout=combineimages(imin) %later with switch? now put side-by-side
samesize=true;
for k=length(imin):-1:1
    s(k,:)=size(imin{k});
    samesize=samesize & all(s(k,:)==s(end,:));
end

if samesize %all same size
    sall=s(1,:);
    sall(2)=sum(s(:,2));
    imout=zeros(sall,'like',imin{1});
    ind=0;
    for k=1:length(imin)
        imout(:,ind+1:ind+s(k,2),:,:)=imin{k};
        ind=ind+s(k,2);
    end
else
    disp('channels of different sizes not yet implemented im imageloaderSMAP')
end
%     imout=imin;
end
