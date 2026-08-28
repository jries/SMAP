classdef imageloaderNDTiff<interfaces.imageloaderSMAP
    %imageloaderNDTiff image loader for NDTiff datasets (Micro-Manager,
    %pycro-manager). An NDTiff dataset is a directory containing
    %NDTiff.index and one or several *_NDTiffStack*.tif files. The images
    %are read directly from the offsets listed in the index, the metadata
    %is parsed from the JSON blocks in the tif files.

    properties
        index           %parsed NDTiff.index, see readNDTiffIndex
        fids            %file identifiers, one per tif file of the dataset
        folder          %directory of the dataset
        summarytext=''  %summary metadata (JSON) of the first tif file
        readoutimgtags={};
        imtags
        init
    end

    methods
        function obj=imageloaderNDTiff(varargin)
            obj@interfaces.imageloaderSMAP(varargin{:});
        end
        function openi(obj,file)
            closefiles(obj);
            if isfolder(file)
                obj.folder=file;
            else
                if ~exist(file,'file')
                    return
                end
                obj.folder=fileparts(file);
            end
            obj.index=readNDTiffIndex(obj.folder);
            obj.file=[obj.folder filesep obj.index.files{obj.index.fileindex(1)}];
            openfiles(obj);
            obj.summarytext=readsummarytext(obj);
            obj.getmetadata;
            obj.metadata.basefile=obj.folder;
            initimagetags(obj);
        end
        function image=getimagei(obj,frame)
            image=[];
            if isempty(obj.index)||isempty(obj.fids)||frame<1||frame>obj.index.n
                return
            end
            fid=obj.fids(obj.index.fileindex(frame));
            if fid<0
                return
            end
            w=obj.index.width(frame);h=obj.index.height(frame);
            if obj.index.bytesPerPixel(frame)==1
                prec='*uint8';
            else
                prec='*uint16';
            end
            if fseek(fid,obj.index.pixelOffset(frame),'bof')~=0
                return
            end
            imh=fread(fid,[w h],prec);
            if numel(imh)<w*h
                image=[];
                return
            end
            image=imh';
            readimagetags(obj,frame);
        end
        function prefit(obj)
            initimagetags(obj);
        end
        function closei(obj)
            closefiles(obj);
        end

        function image=getimageonline(obj,number)
            image=obj.getimage(number);
            if isempty(image)&&obj.onlineAnalysis
                disp('wait')
                pause(obj.waittime*2)
                reloadindex(obj);
                image=obj.getimage(number);
            end
        end

        function allmd=getmetadatatagsi(obj)
            allmd=parsejsonflat(readmetadatatext(obj,1)); %image metadata first, as in imageloaderMM
            allmd=vertcat(allmd,{'Format','NDTiff'});

            %ROI, same format as in the micro manager metadata: 'x-y-w-h'
            try
                roistr=gettagi(allmd,'ROI');
                if ~isempty(roistr)
                    troi=textscan(roistr,'%d','delimiter','-');
                    allmd(end+1,:)={'ROI direct',num2str(troi{:}')};
                end
            catch err
                disp(getReport(err,'basic'))
            end

            allmd(end+1,:)={'frames direct',num2str(obj.index.n)};
            allmd(end+1,:)={'Width info',obj.index.width(1)};
            allmd(end+1,:)={'Height info',obj.index.height(1)};

            %time between frames from the elapsed time of the first and last image
            dt=NaN;
            try
                if obj.index.n>1
                    mdlast=parsejsonflat(readmetadatatext(obj,obj.index.n));
                    t0=str2double(gettagi(allmd,'ElapsedTime-ms'));
                    te=str2double(gettagi(mdlast,'ElapsedTime-ms'));
                    dt=(te-t0)/(obj.index.n-1);
                end
            catch err
                disp(getReport(err,'basic'))
            end
            allmd(end+1,:)={'timediff direct',num2str(dt)};

            try
                txt=fileread([obj.folder filesep 'comments.txt']);
                cfile=jsondecode(txt);
                ctxt=cfile.map.GeneralAnnotation.scalar.comments.scalar;
                allmd(end+1,:)={'Comments direct',strrep(ctxt,newline,'; ')};
            catch
            end

            allmd=vertcat(allmd,parsejsonflat(obj.summarytext));
            obj.allmetadatatags=allmd;
        end
    end
end

function openfiles(obj)
obj.fids=zeros(length(obj.index.files),1);
for k=1:length(obj.index.files)
    obj.fids(k)=fopen([obj.folder filesep obj.index.files{k}],'r','l');
end
end

function closefiles(obj)
for k=1:length(obj.fids)
    if obj.fids(k)>0
        fclose(obj.fids(k));
    end
end
obj.fids=[];
end

function reloadindex(obj)
%during online analysis the index and the tif files keep growing
try
    closefiles(obj);
    obj.index=readNDTiffIndex(obj.folder);
    openfiles(obj);
    obj.metadata.numberOfFrames=obj.index.n;
catch err
    disp(getReport(err,'basic'))
end
end

function txt=readsummarytext(obj)
%the summary metadata is stored in the header of the tif file:
%bytes 20-23: 2355492 (marker), bytes 24-27: length, then the JSON
txt='';
fid=obj.fids(1);
if fid<0
    return
end
try
    fseek(fid,20,'bof');
    marker=fread(fid,1,'*uint32');
    len=double(fread(fid,1,'*uint32'));
    if marker==2355492&&len>0
        txt=fread(fid,len,'*char')';
    end
catch err
    disp(getReport(err,'basic'))
end
end

function txt=readmetadatatext(obj,frame)
txt='';
if isempty(obj.index)||isempty(obj.fids)||frame<1||frame>obj.index.n
    return
end
fid=obj.fids(obj.index.fileindex(frame));
if fid<0
    return
end
if fseek(fid,obj.index.mdOffset(frame),'bof')~=0
    return
end
txt=fread(fid,obj.index.mdLength(frame),'*char')';
txt=deblank(strrep(txt,char(0),' '));
end

function initimagetags(obj)
camset=obj.getPar('loc_cameraSettings');
if ~isempty(camset)&& myisfield(camset,'imagemetadata') && ~isempty(camset.imagemetadata)
    if iscell(camset.imagemetadata)
        obj.readoutimgtags=camset.imagemetadata;
    else
        obj.readoutimgtags=split(camset.imagemetadata,',');
    end
    obj.readoutimgtags=strtrim(obj.readoutimgtags);
    obj.init=true;
else
    obj.readoutimgtags={};
end
end

function readimagetags(obj,frame)
if isempty(obj.readoutimgtags)
    return
end
if obj.init
    obj.imtags=zeros(length(obj.readoutimgtags),obj.metadata.numberOfFrames);
    obj.init=false;
end
txt=readmetadatatext(obj,frame);
for k=1:length(obj.readoutimgtags)
    obj.imtags(k,frame)=str2double(getjsonvalue(txt,obj.readoutimgtags{k}));
end
end

function val=getjsonvalue(txt,key)
%read out a single key of a flat JSON string without parsing all of it
val='';
if isempty(txt)||isempty(key)
    return
end
n=length(txt);
ind=strfind(txt,['"' key '"']);
for m=1:length(ind)
    p=ind(m)+length(key)+2;
    while p<=n&&isspace(txt(p))
        p=p+1;
    end
    if p>n||txt(p)~=':' %key of a different tag or a value, not a name
        continue
    end
    p=p+1;
    while p<=n&&isspace(txt(p))
        p=p+1;
    end
    if p>n
        return
    end
    if txt(p)=='"'
        [val,~]=readjsonstring(txt,p);
    else
        pe=p;
        while pe<=n&&txt(pe)~=','&&txt(pe)~='}'
            pe=pe+1;
        end
        val=strtrim(txt(p:pe-1));
    end
    return
end
end

function val=gettagi(md,tag)
val='';
if isempty(md)
    return
end
ind=find(strcmp(md(:,1),tag),1,'first');
if ~isempty(ind)
    val=md{ind,2};
end
end

function md=parsejsonflat(txt)
%parse the top level of a JSON object into a n x 2 cell array of names and
%values. Values are returned as strings (as gethashtable does for the
%java metadata of imageloaderMM), nested objects and arrays are skipped.
%The key names are preserved as they are (jsondecode would mangle names
%like 'Andor-Camera', which are used to identify the camera).
md=cell(0,2);
if isempty(txt)
    return
end
n=length(txt);
k=find(txt=='{',1,'first');
if isempty(k)
    return
end
k=k+1;
mdh=cell(1000,2);
nmd=0;
while k<=n
    while k<=n&&(isspace(txt(k))||txt(k)==',')
        k=k+1;
    end
    if k>n||txt(k)=='}'
        break
    end
    if txt(k)~='"' %unexpected, skip
        k=k+1;
        continue
    end
    [key,k]=readjsonstring(txt,k);
    while k<=n&&(isspace(txt(k))||txt(k)==':')
        k=k+1;
    end
    if k>n
        break
    end
    switch txt(k)
        case '"'
            [val,k]=readjsonstring(txt,k);
        case {'{','['}
            k=skipjsonvalue(txt,k);
            continue
        otherwise
            ks=k;
            while k<=n&&txt(k)~=','&&txt(k)~='}'
                k=k+1;
            end
            val=strtrim(txt(ks:k-1));
    end
    nmd=nmd+1;
    if nmd>size(mdh,1)
        mdh(2*size(mdh,1),2)={[]};
    end
    mdh(nmd,1:2)={key,val};
end
md=mdh(1:nmd,:);
end

function [str,k]=readjsonstring(txt,k)
%txt(k) is the opening ", returns the unescaped string and the position
%after the closing "
n=length(txt);
k=k+1;
ks=k;
esc=false;
while k<=n
    c=txt(k);
    if c=='\'
        esc=true;
        k=k+2;
        continue
    elseif c=='"'
        break
    end
    k=k+1;
end
str=txt(ks:min(k,n+1)-1);
k=k+1;
if esc
    str=strrep(str,'\n',newline);
    str=strrep(str,'\r',char(13));
    str=strrep(str,'\t',char(9));
    str=strrep(str,'\/','/');
    str=strrep(str,'\"','"');
    str=strrep(str,'\\','\');
end
end

function k=skipjsonvalue(txt,k)
%skip a nested object or array starting at txt(k)
n=length(txt);
depth=0;
while k<=n
    c=txt(k);
    if c=='"'
        [~,k]=readjsonstring(txt,k);
        continue
    elseif c=='{'||c=='['
        depth=depth+1;
    elseif c=='}'||c==']'
        depth=depth-1;
        if depth<=0
            k=k+1;
            return
        end
    end
    k=k+1;
end
end
