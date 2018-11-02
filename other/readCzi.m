function [out,meta]=readCzi(imin)
data=imin{1};
ss=size(data{1,1});
out=zeros(ss(1),ss(2),length(data)/2,2);
for k=1:length(data)/2
    out(:,:,k,1)=data{2*(k-1)+1,1};
    out(:,:,k,2)=data{2*k,1};
end
md=imin{2};
% md.get('Global Information|Image|Channel|LaserScanInfo|PixelTime #1');
meta.frametime=str2double(md.get('Global Information|Image|Channel|LaserScanInfo|FrameTime #1'));
meta.totaltime=meta.frametime*length(data)/2;