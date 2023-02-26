global pn
framerate=15;

[fn,pn]=uigetfile([pn filesep '*.*']);

[~,~,ext]=fileparts(fn);
switch ext
    case '.tif'
        video=tiffreadVolume([pn filesep fn]);
        video=permute(video,[1,2,4,3]);
    case '.gif'
        [img,map]=imread([pf f],'Frames','all');
        img=squeeze(img)*intensityscaling;
        clear video;
        for k=1:size(img,3)
            video(k)=im2frame(img(:,:,k),map);
        end
    otherwise
        v=VideoReader([pn filesep fn]);
        video=read(v);
        framerate=v.FrameRate;
end
% %convert 2 gray
% saturation=1.5;
% video=(sum(video,3));
% video=uint8(video/max(video(:))*saturation*2^8);

fout=[pn filesep strrep(fn,ext,'_c.mp4')];
mysavemovie(video,fout,'FrameRate',framerate)