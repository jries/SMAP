%convert2mp4
framerate=5;

%%
global pf
[f,pf]=uigetfile([pf '*.*']);
[~,fn,ext]=fileparts(f);

%%
switch ext
    case '.gif'
        [img,map]=imread([pf f],'Frames','all');
        img=squeeze(img);
    otherwise
        disp('not yet implemented')
end

clear fr;
for k=1:size(img,3)
    fr(k)=im2frame(img(:,:,k),map);
end

%%
mysavemovie(fr,[pf fn '_n.mp4'],'FrameRate',framerate)



