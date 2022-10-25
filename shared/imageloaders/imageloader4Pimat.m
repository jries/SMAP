classdef imageloader4Pimat<interfaces.imageloaderSMAP
    %imageloader4Pimat image loader for 4Pi image files, after conversion
    %to mat
    properties
         reader
         batchlength
         allfiles
         currentfiledata %.q1,..., .batch
    end
    
    methods
        function obj=imageloader4Pimat(varargin)
            obj@interfaces.imageloaderSMAP(varargin{:});
        end
        function openi(obj,file)
            % find all files in folder
            path=fileparts(file);
            allfiles=dir([path filesep '*.mat']);
            obj.allfiles.files={allfiles(:).name};
            obj.allfiles.path=path;
            md=obj.getmetadata;
            getmetadatatagsi(obj);
            obj.metadata.basefile=path;
            obj.currentfiledata.batch=0;
            obj.currentfiledata.frames=[];
            obj.currentImageNumber=1;
            obj.file=file;
        end
        function image=getimagei(obj,frame)
%             try
                batch=ceil(frame/obj.batchlength);
                fnum=mod(frame-1,obj.batchlength)+1;
                if batch~=obj.currentfiledata.batch
                    if batch>length(obj.allfiles.files) 
                        image=[];
                        return;
                    end
                    obj.currentfiledata.frames=load([obj.allfiles.path filesep obj.allfiles.files{batch}]);
                    obj.currentfiledata.batch=batch;
                end
                if fnum<=size(obj.currentfiledata.frames.qd1,3)
                    image(:,:,1,1)=obj.currentfiledata.frames.qd1(:,:,fnum);
                    image(:,:,1,2)=obj.currentfiledata.frames.qd2(:,:,fnum);
                    image(:,:,1,3)=obj.currentfiledata.frames.qd3(:,:,fnum);
                    image(:,:,1,4)=obj.currentfiledata.frames.qd4(:,:,fnum);
                else
                    image=[];
                end

%                 image={obj.currentfiledata.frames.qd1(:,:,fnum),obj.currentfiledata.frames.qd2(:,:,fnum),obj.currentfiledata.frames.qd3(:,:,fnum),obj.currentfiledata.frames.qd4(:,:,fnum)};
%             catch err
%                 err
%                 image=[];
%             end
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
        
        function allmd=getmetadatatagsi(obj)
            allmd={'Format','4PiMat'};
                        % determine batch length
            f1=[obj.allfiles.path filesep obj.allfiles.files{1}];
            f1info=whos('-file',f1);
            allmd(end+1,:)={'batchlength',f1info(1).size(3)};
            allmd(end+1,:)={'NumberOfChannels',4};
            obj.batchlength=f1info(1).size(3);


            allmd(end+1,:)={'Width info',f1info(1).size(2)};
            allmd(end+1,:)={'Height info',f1info(1).size(1)};

            allmd(end+1,:)={'FileName',f1};
            allmd(end+1,:)={'Roi direct',num2str([0 0 f1info(1).size(2) f1info(1).size(1)])};

            fe=[obj.allfiles.path filesep obj.allfiles.files{end}];
            feinfo=whos('-file',fe);
            numf=(length(obj.allfiles.files)-1)*obj.batchlength+feinfo(1).size(3);
            allmd(end+1,:)={'Frames',numf};
            allmd(end+1,:)={'frames direct',numf};
            obj.allmetadatatags=allmd;
        end       
    end    
end