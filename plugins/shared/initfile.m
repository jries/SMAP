function filestruct=initfile(filename)
tifs=struct('image',[],'info',[]);
% finfo=obj.getPar('loc_fileinfo');
% sfile=finfo.basefile;
if nargin<1||isempty(filename)
    filename='smapgenerated';
end
if isempty(strfind(filename,'_sml'))
    filename=[filename '_sml.mat'];  
end
infost=struct('camId','SMAP generated file','port','none','exposure',0,'emgain',1,...
    'conversion',1,'offset',0,'cam_pixelsize_um',[.1 .1],'roi',[0  0 512 512], 'comment','',...
    'filename',filename,'numberOfFrames',0,'Width',512,'Height',512,'format','none');


%git version
[gitstatus,head]=system('git rev-parse HEAD');
if gitstatus == 0
infost.git.head=head;
else
    infost.git.head=[];
end
[gitstatus,branch]=system('git status');
if gitstatus== 0
ind=strfind(branch,'On branch');
le=find(branch(ind+1:end)==10,1,'first')+ind;
infost.git.branch=branch(ind+10:le-1);
else
infost.git.branch='not found';
end

filestruct=struct('info',infost,'average',[],'name',filename,...
    'number',1,...
    'numberOfTif',0,'tif',tifs,'raw',struct('image',[],'frame',[]));
end
