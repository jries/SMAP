function mysavemovie(img,fsave,varargin)
% varargin: name, value pairs to be passed on to movie constructor
% img: frames or array
% fsave: output path

v=VideoWriter([pf fout],'MPEG-4');
v.FrameRate=obj.getSingleGuiParameter('framerate');
open(v);
writeVideo(v,obj.outmovie(1:end-1))
close(v);