classdef imageloader_mat<interfaces.imageloaderSMAP
    %imageloader4Pimat image loader for 4Pi image files, after conversion
    %to mat
    properties
         reader
         batchlength
         infofiles
         infoall
         allfiles
         currentfiledata %.q1,..., .batch
         lastimagepos
    end
    
    methods
        function obj=imageloader_mat(varargin)
            obj@interfaces.imageloaderSMAP(varargin{:});
        end
        function openi(obj,file)
            % find all files in folder
            file=char(file);
            %extract pattern:
            ind=strfind(file,'_');
            path=fileparts(file);
            searchstr=[file(1:ind(end-1)) '*' file(ind(end):end)];
            allf=dir(searchstr);
            obj.allfiles.files={allf(:).name};
            obj.allfiles.path=path;
            analyzecontent(obj);
            md=obj.getmetadata;
            getmetadatatagsi(obj);
            obj.metadata.basefile=path;
            obj.currentfiledata.batch=0;
            obj.currentfiledata.frames=[];
            obj.currentImageNumber=1;
            obj.file=file;
        end
        function image=getimagei(obj,frame)
            batch=find(frame<=obj.infoall.framesinbatch,1,'last');
            if batch>1
                fnum=frame-obj.infoall.framesinbatch(batch-1);
            else
                fnum=frame;
            end
            if batch~=obj.currentfiledata.batch
                if batch>length(obj.allfiles.files) 
                    image=[];
                    return;
                end
                obj.currentfiledata.frames=load([obj.allfiles.path filesep obj.allfiles.files{batch}]);
                obj.currentfiledata.batch=batch;
            end
            image=[];
            [tind,zind]=ind2sub([obj.infoall.timepoints,obj.infoall.zslices],fnum);
            if tind>obj.infoall.timepoints || zind>obj.infoall.zslices
                return
            end

            for k=obj.infoall.channels:-1:1
                    image(:,:,1,k)=obj.currentfiledata.frames.(obj.infofiles(batch).name{k})(:,:,tind,zind);
            end
            obj.lastimagepos=[batch, tind, zind];
        end
        
        function closei(obj)
        end
        
        function image=getimageonline(obj,number)
            image=obj.getimage(number);
            if isempty(image)&&obj.onlineAnalysis 
                    disp('wait')
                    pause(obj.waittime*2)
                    image=obj.getimage(number);
            end
        end
        function analyzecontent(obj)
            numframes=0;
            for k=1:length(obj.allfiles.files)
                f1=[obj.allfiles.path filesep obj.allfiles.files{k}];
                f1info=whos('-file',f1);
                info(k).channels=0;
                indch=1;
                for l=1:length(f1info)
                    if strcmpi(f1info(l).name,'metadata')
                    elseif contains(f1info(l).name,'qd')
                        sh=f1info(l).size;
                        info(k).channels=info(k).channels+1;
                        info(k).zslices(indch)=sh(4);
                        info(k).timepoints(indch)=sh(3);
                        info(k).name{indch}=f1info(l).name;
                        indch=indch+1;
                    else
                        disp(f1info(l).name)
                    end
                end
                numframes=numframes+max(info(k).zslices)*max(info(k).timepoints);
                zslicesall(k)=max(info(k).zslices); timepointsall(k)=max(info(k).timepoints); 
            end
            obj.infoall.numframes=numframes;
            obj.infoall.zslices=zslicesall;
            obj.infoall.timepoints=timepointsall;
            obj.infoall.framesinbatch=cumsum(zslicesall*timepointsall);
            obj.infoall.channels=max([info(:).channels]);
            obj.infofiles=info;
            obj.infoall.imagesize=sh(1:2);
        end

        function allmd=getmetadatatagsi(obj)
            allmd={'Format','Mat'}; 
            allmd(end+1,:)={'NumberOfChannels',obj.infoall.channels};
            allmd(end+1,:)={'ZSlices',obj.infoall.zslices};
            allmd(end+1,:)={'Timepoints',obj.infoall.timepoints};
            allmd(end+1,:)={'Width info',obj.infoall.imagesize(2)};
            allmd(end+1,:)={'Height info',obj.infoall.imagesize(1)};
            f1=[obj.allfiles.path filesep obj.allfiles.files{1}];
            allmd(end+1,:)={'FileName',f1};
            allmd(end+1,:)={'Roi direct',num2str([0 0 obj.infoall.imagesize(2) obj.infoall.imagesize(1)])};
            allmd(end+1,:)={'Frames',obj.infoall.numframes};
            allmd(end+1,:)={'frames direct',obj.infoall.numframes};
            

            l=load(f1,'metadata');
            if isfield(l,'metadata')
                list=mystruct2list(l.metadata);
                allmd=vertcat(allmd,list);
            end

            obj.allmetadatatags=allmd;
        end       
    end    
    end

