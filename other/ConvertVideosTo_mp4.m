global pn

[fn,pn]=uigetfile([pn filesep '*.*']);

[~,~,ext]=fileparts(fn);
switch ext
    case '.tif'

        video=tiffreadVolume([pn filesep fn]);
        video=permute(video,[1,2,4,3]);
    otherwise
        v=VideoReader([pn filesep fn]);
        video=read(v);
end
% %convert 2 gray
% saturation=1.5;
% video=(sum(video,3));
% video=uint8(video/max(video(:))*saturation*2^8);


fout=[pn filesep strrep(fn,ext,'_c.mp4')];
mysavemovie(video,fout,'FrameRate',v.FrameRate)