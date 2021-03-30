global pn

[fn,pn]=uigetfile([pn filesep '*.avi']);

v=VideoReader([pn filesep fn]);

video=read(v);
fout=[pn filesep strrep(fn,'.avi','_c.mp4')];
mysavemovie(video,fout,'FrameRate',v.FrameRate)