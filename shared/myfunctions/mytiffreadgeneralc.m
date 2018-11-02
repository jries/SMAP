function [imag,info,fileout]=mytiffreadgeneralc(file,frames,opt,calfile)
if nargin<4||isempty(calfile)
        calfile='CameraCalibration.xls';
end
%info: gain, size, all...
[pathstr, name, ext]=fileparts(file);

if isempty(name)
    df=myfastdir(pathstr,  'img*.tif');
    whos df
    name=df(1);
    file=[pathstr filesep name];
end

    dind=findstr(name,'_');
searchx=[name(1:dind(1)) '*' ext];
allfiles=myfastdir(pathstr,searchx);
    fileout.names=allfiles;
    if frames==0
        frames=1:length(allfiles);
    end
    fileout.pfad=[pathstr filesep];
    fileout.frames=frames;
    fileout.numframes=length(allfiles);
    

    finfo=imfinfo(file);
    [p n x]=fileparts(file);
    
    
    metafile=([ p filesep 'metadata.txt'])
    fid=fopen(metafile,'r');
    if fid>0 %MM metadata present
        fclose(fid);
        
        camcalib=readtable(calfile);
        minfo=fileread(metafile);
        info=minfoparsec(minfo,camcalib);

    else %no metafile present: now: assign dummy values. later read sif etc.
        mf=dir([pathstr filesep name(1:dind(1)-1) '*.sif'])
        mff=[pathstr filesep mf.name];
        fid=fopen(mff);
        if fid>0 %Andor SIF file
            fclose(fid);
            info=mysifinfo2(mff);
        else
        disp('no metafile present. Default values.')
%             info.conversion=12;
%             info.actualconv=10.9;
%             info.conversion=10.88;
%             info.actualconv=10.88;
            info.actualconv=10;
            info.conversion=10;

            info.readoutrate=1;
            info.exposure=1;
            info.gain=1;
            info.offset=0;
            info.emgain=100;  
            info.chip='not determined';
            info.timediff=1;
            info.port='EM';
            info.roi=[];
            info.cam_pixelsize_um=0.1;  %depends on settings, camera
        end

    end
    


info.Width= finfo(1).Width;
info.Height=finfo(1).Height;


w=finfo(1).Width;
h=finfo(1).Height;
finfo(1).BitDepth
if finfo(1).BitDepth==16
    vtype='uint16';
elseif finfo(1).BitDepth==8
    vtype='uint8';
end

if nargin>=3
    if strcmp(opt,'numf')
        dx=ceil(fileout.numframes/frames);
        frames=ceil(dx/2):dx:fileout.numframes;
    end
end

nframes=length(frames);
entries=h*w*nframes;

% [user,sys] = memory;
% free=sys.PhysicalMemory.Available
free=3e9;
maxmem=free*0.7;
if entries*4>maxmem %4GB max
    nframes=round(maxmem/h/w/4);
    errordlg(['too much memory required, loading only first ' num2str(frames) ' frames'])

end

imag=zeros(h,w,nframes,vtype); 

fmt_s = imformats('tif');

for k=1:nframes
    if ~abs((round(k/10)-k/10))
       
    end
    if frames(k)<=length(allfiles)
    fr=myimread([pathstr filesep allfiles{frames(k)}],fmt_s);
%         imag(:,:,k)=single(fr);
        imag(:,:,k)=(fr);
        
%         figure(1)
%         imagesc(fr)
%         asdf
    else
        break;
    end
end

info.frames=frames;
fileout.frames=frames;
info.size=[h w frames];

% if ~mm
%     info.offset=min(imag(:));
% end
