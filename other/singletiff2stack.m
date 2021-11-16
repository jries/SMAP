%singletiff2stack
global path
[fn, path]=uigetfile([path filesep '*.tif'],'MultiSelect','on');

[~,sind]=sort(fn);

img=imread([path fn{1}]);

imstack=zeros(size(img,1),size(img,2),length(sind),'like',img);

for k=1:length(sind)
    imstack(:,:,k)=imread([path fn{sind(k)}]);
end

[~,imgname]=fileparts(fn{1});
outname=[path(1:end-1) '_' imgname '_combined.tif'];

saveastiff(imstack,outname)
