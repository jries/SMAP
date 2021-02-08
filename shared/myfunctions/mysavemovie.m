function mysavemovie(img,fsave,varargin)
% varargin: name, value pairs to be passed on to movie constructor
% look at documentation from VideoWriter
% 'profile': standard VideoWriter profile. Default: MPEG-4
% img: frames or array
% fsave: output path

% framerate
ind=find(strcmp(varargin,'profile'),1,'first');
if ~isempty(ind)
    profile=varargin{ind+1};
    varargin(ind:ind+1)=[];
else
    profile='MPEG-4';
end

v=VideoWriter(fsave,profile);

for k=1:2:length(varargin)
    set(v,varargin{k},varargin{k+1});
end

open(v);
writeVideo(v,img)
close(v);