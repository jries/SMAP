function [imstack, roi, pixelsize,settings3D]=readbeadimages(file,p)
[path,f,ext]=fileparts(file);
indq=strfind(f,'_q');
settings3D=[];
multichannel_4pi=false;
if ~isempty(indq)
    allfiles=dir([path filesep f(1:indq+1)  '*' ext]);
    for k=1:length(allfiles)
        files{k}=[path filesep allfiles(k).name];
        disp(allfiles(k).name(indq:end))
    end
    file=files;
    multichannel_4pi=true;
end
pixelsize=100;
    
    if isfield(p,'smap') && p.smap
        
        try
             r=imageloaderAll(file,[],p.smappos.P);

             imstack=r.getmanyimages(1:r.metadata.numberOfFrames,'mat');
             if size(imstack,3)<=1 %only one frame read
                 error('only one frame found')
             end
             imstack=imstack-r.metadata.offset;
             roi=r.metadata.roi;
             pixelsize=r.metadata.cam_pixelsize_um;
             r.close;
        catch err
            err
            imstack=readfile_tif(file);
            roi=[0 0 size(imstack,1) size(imstack,2)]; %check x,y
        end
        if isempty(imstack)
            disp('using simple reader')
            warndlg('using simple reader, this might create problems if only part of the camera chip is used.','using simple reader','replace');
            if multichannel_4pi
                imstack=[];
                for k=1:length(file)
                    imstack=horzcat(imstack,readfile_tif(file{k}));
                end
            else
                 imstack=readfile_tif(file);
            end
            roi=[0 0 size(imstack,1) size(imstack,2)];
        end
          
    else
        imstack=readfile_tif(file);
        roi=[0 0 size(imstack,1) size(imstack,2)]; %check x,y
    end
    if multichannel_4pi
        wx=size(imstack,2)/4;wy=size(imstack,1);
        settings3D=struct('y4pi',[0 0 0 0],'x4pi',[0 wx 2*wx 3*wx], 'width4pi',wx,'height4pi',wy,'mirror4pi',[0 0 0 0],'pixelsize_nm',100,'offset',100,'conversion',0.5);
    end
    if isfield(p,'framerange')
        if length(p.framerange)~=2
            fr=p.framerange;
        else
            fr=p.framerange(1):p.framerange(2);
        end
           
        imstack=imstack(:,:,fr);
    end
    if isfield(p,'emgain') && p.emgain
        imstack=imstack(:,end:-1:1,:);
        %XXXX for Andor take out?
%          if any(roi(1:2)>0) %if roi(1:2)=[0 0] it is likely that roi was not read out and set to default.
%              roi(1)=512-roi(1)-roi(3);
%          end
    end
%     if isfield(p,'roimask')&&~isempty(p.roimask)
%         imstack=imstack.*uint16(p.roimask);
%     end
end