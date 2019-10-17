classdef imageloaderMMsingle<interfaces.imageloaderSMAP
    %imageloaderMM image loader for micromanager single tiff files
    %   Detailed explanation goes here
    
    properties
%         calfile='settings/CameraCalibration.xls';
        separate
        separatefiles
    end
    
    methods
        function obj=imageloaderMMsingle(varargin)
            obj@interfaces.imageloaderSMAP(varargin{:});
        end
        function openi(obj,file)
            obj.file=file;
            info=getimageinfo(file);
            obj.getmetadata;
            obj.metadata.basefile=info.basefile;
            
            obj.separatefiles=info.files;
            obj.separate.numfiles=info.frames;
            obj.separate.path=info.path;
%             obj.currentImageNumber=0;
            obj.separate.fmt_s=imformats('tif');
        end
%         function mdo=getmetadata(obj)
%             mdo=getmetadata@interfaces.imageloaderSMAP(obj);
%             obj.metadata.basefile=info.basefile;
% %             mdo2=getmetadataMM(obj); 
%             
%             
%         end
        function allmd=getmetadatatagsi(obj)
%             metafile=[fileparts(obj.file) filesep 'metadata.txt'];
            allmd=getmetadataMMnew(obj.file);
        end
        function image=getimagei(obj,frame)
            image=readseparate(obj,frame);
        end
        function closei(obj)
            %not needed, as single files are read
        end
        function file= getbasefile(obj)
            file=getbasefile(obj.file);
        end
    end
    
end

function image=readseparate(obj,number)
image=[];
separate=obj.separate;
% lenfiles=length(separate.files);
lenfiles=separate.numfiles;
if lenfiles<number % || (ismethod(obj.separatefiles(number),'ismissing') && obj.separatefiles(number).ismissing)
    if obj.onlineAnalysis %ask for image that is not in list
        lastfile= obj.separatefiles{lenfiles};
%         thisname= generateFileName(lastfile,lenfiles,obj.metadata.allmetadata.numberNameRange,number);
        thisname= generateFileName(lastfile,lenfiles,number);
        thisfile=[separate.path filesep thisname];
        if ~exist(thisfile,'file')
            disp('wait')
            pause(obj.waittime*2)
        end
        if ~exist(thisfile,'file')
            image=[];
            return
        else
%             if number>lenfiles
%                 obj.separatefiles{number+1000}='';
%             end
            obj.separatefiles{number}=thisname;
            obj.separate.numfiles=max(lenfiles,number);
            obj.metadata.numberOfFrames=max(lenfiles,number);
        end 
    else
        image=[];
        return
    end
else
    try
        
    thisfile=[separate.path filesep obj.separatefiles{number}];
    catch err
        
        err
        return
    end
end
try
    image=myimread(thisfile,separate.fmt_s);     
catch err
    pause(1)
    disp('error encountered rading image, try again');
    image=myimread(thisfile,separate.fmt_s);     
end
end

function newfile=generateFileName(oldfile, oldfilenumber,newfilenumber)
% if isempty(indfbar)
    oldfn=[num2str(oldfilenumber-1) '_'];
    inds=strfind(oldfile,oldfn);
    if strcmp(oldfile(inds-1),'0')
        oldfn=['0' oldfn];
    end
    newfn=[num2str(newfilenumber-1) '_'];
    if length(newfn)<length(oldfn)
        newfn=['0' newfn];
    end
    newfile=strrep(oldfile,oldfn,newfn);
% end
% oldfilenamenumber=str2double(oldfile(indfbar(1):indfbar(2)));
% newfilenamenumber=oldfilenamenumber-oldfilenumber+newfilenumber;
% lenfield=indfbar(2)-indfbar(1)+1;
% newfilestr=num2str(newfilenamenumber,['%0' num2str(lenfield) 'i']);
% newfile=[oldfile(1:indfbar(1)-1) newfilestr oldfile(indfbar(2)+1:end)];
end