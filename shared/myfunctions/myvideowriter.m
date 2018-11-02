function myvideowriter(mov,filename)
% video writer

myVideo=VideoWriter(filename);
open(myVideo)
writeVideo(myVideo,mov)
close(myVideo);