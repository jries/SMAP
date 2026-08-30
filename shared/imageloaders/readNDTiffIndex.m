function index=readNDTiffIndex(folder)
% readNDTiffIndex parses the NDTiff.index file of an NDTiff dataset
% (Micro-Manager / pycro-manager NDTiffStorage).
% The index is a flat binary table with one variable-length record per
% image, all fields little-endian uint32:
%   axesLength, axes (JSON), filenameLength, filename,
%   pixelOffset, width, height, pixelType, pixelCompression,
%   metadataOffset, metadataLength, metadataCompression
% The table is zero-padded at the end, parsing stops at the first record
% with axesLength==0.
% Returns a struct with one array per field, sorted by the time axis.

indexfile=[folder filesep 'NDTiff.index'];
if ~exist(indexfile,'file')
    error('readNDTiffIndex:noindex','no NDTiff.index found in %s',folder);
end

fid=fopen(indexfile,'r','l');
d=fread(fid,inf,'*uint8');
fclose(fid);
n=length(d);

maxentries=ceil(n/50)+16;
pixelOffset=zeros(maxentries,1);
width=zeros(maxentries,1);
height=zeros(maxentries,1);
pixelType=zeros(maxentries,1);
compression=zeros(maxentries,1);
mdOffset=zeros(maxentries,1);
mdLength=zeros(maxentries,1);
fileindex=zeros(maxentries,1);
time=zeros(maxentries,1);
files={};
lastname='';lastind=0;
hasextraaxes=false;

p=1;k=0;
while p+3<=n
    alen=double(typecast(d(p:p+3),'uint32'));p=p+4;
    if alen==0 %zero padding at the end of the index
        break
    end
    if p+alen-1>n
        break
    end
    ax=char(d(p:p+alen-1).');p=p+alen;
    if p+3>n
        break
    end
    flen=double(typecast(d(p:p+3),'uint32'));p=p+4;
    if flen==0||p+flen-1>n
        break
    end
    fname=char(d(p:p+flen-1).');p=p+flen;
    if p+31>n
        break
    end
    v=double(typecast(d(p:p+31),'uint32'));p=p+32;

    k=k+1;
    if k>maxentries %should not happen, but do not fail on it
        maxentries=2*maxentries;
        pixelOffset(maxentries)=0;width(maxentries)=0;height(maxentries)=0;
        pixelType(maxentries)=0;compression(maxentries)=0;mdOffset(maxentries)=0;
        mdLength(maxentries)=0;fileindex(maxentries)=0;time(maxentries)=0;
    end
    pixelOffset(k)=v(1);width(k)=v(2);height(k)=v(3);pixelType(k)=v(4);
    compression(k)=v(5);mdOffset(k)=v(6);mdLength(k)=v(7);

    if ~strcmp(fname,lastname) %usually all entries point to the same file
        lastind=find(strcmp(files,fname),1,'first');
        if isempty(lastind)
            files{end+1}=fname; %#ok<AGROW>
            lastind=length(files);
        end
        lastname=fname;
    end
    fileindex(k)=lastind;

    it=strfind(ax,'"time":');
    if isempty(it)
        time(k)=NaN;
    else
        th=sscanf(ax(it(1)+7:end),'%d',1);
        if isempty(th)
            time(k)=NaN;
        else
            time(k)=th;
        end
    end
    if length(strfind(ax,':'))>1
        hasextraaxes=true;
    end
end

if k==0
    error('readNDTiffIndex:empty','no valid entries found in %s',indexfile);
end

ind=1:k;
pixelOffset=pixelOffset(ind);width=width(ind);height=height(ind);
pixelType=pixelType(ind);compression=compression(ind);mdOffset=mdOffset(ind);
mdLength=mdLength(ind);fileindex=fileindex(ind);time=time(ind);

if any(compression~=0)
    error('readNDTiffIndex:compressed','compressed NDTiff images are not supported');
end

if hasextraaxes
    disp('NDTiff: dataset has axes other than time. Images are loaded in the order of the index.')
end

if all(~isnan(time))&&~issorted(time)
    [time,sortind]=sort(time);
    pixelOffset=pixelOffset(sortind);width=width(sortind);height=height(sortind);
    pixelType=pixelType(sortind);
    mdOffset=mdOffset(sortind);mdLength=mdLength(sortind);fileindex=fileindex(sortind);
end

bytesPerPixel=2*ones(k,1);
bytesPerPixel(pixelType==0)=1;
%sanity check with the offset of the metadata, which directly follows the
%pixels
dm=mdOffset(1)-pixelOffset(1);
if dm==width(1)*height(1)
    bytesPerPixel(:)=1;
elseif dm==2*width(1)*height(1)
    bytesPerPixel(:)=2;
end

%remove entries that point beyond the end of their tif file (acquisition
%truncated, e.g. when the 4 GB limit of the 32 bit offsets was reached)
filebytes=zeros(length(files),1);
for f=1:length(files)
    dh=dir([folder filesep files{f}]);
    if isempty(dh)
        filebytes(f)=0;
    else
        filebytes(f)=dh.bytes;
    end
end
valid=pixelOffset+width.*height.*bytesPerPixel<=filebytes(fileindex);
if ~all(valid)
    kv=find(~valid,1,'first')-1;
    disp(['NDTiff: index lists ' num2str(k) ' images, but only ' num2str(kv) ...
        ' fit into the tif file(s). Acquisition truncated, using ' num2str(kv) ' frames.'])
    ind=1:kv;
    pixelOffset=pixelOffset(ind);width=width(ind);height=height(ind);
    pixelType=pixelType(ind);bytesPerPixel=bytesPerPixel(ind);
    mdOffset=mdOffset(ind);mdLength=mdLength(ind);fileindex=fileindex(ind);
    time=time(ind);
    k=kv;
end

index.n=k;
index.files=files;
index.fileindex=fileindex;
index.pixelOffset=pixelOffset;
index.width=width;
index.height=height;
index.pixelType=pixelType;
index.bytesPerPixel=bytesPerPixel;
index.mdOffset=mdOffset;
index.mdLength=mdLength;
index.time=time;
index.folder=folder;
end
