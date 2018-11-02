function [b,p]=images2beads(obj,p)
sf=selectManyFiles(p.pathhere);
pos=sf.guihandles.freepos.Position;
sf.guihandles.filtersize=uicontrol('Style','edit','String','3','Parent',sf.handle,'Position',pos);
pos(2)=pos(2)+30;
uicontrol('Style','text','String','filtersize:','Parent',sf.handle,'Position',pos);

waitfor(sf.handle);
filelist=sf.filelist;
b=[];
fs=str2double(sf.addpar.filtersize.String);
h=fspecial('gaussian',round(fs*3),fs);
fmax=0;
roisizeh=round(1.5*(p.roisize-1)/2);
rsr=-roisizeh:roisizeh;
for k=1:length(filelist)
    il=getimageloader(obj,filelist{k});
    imstack=il.getmanyimages([],'mat');
    mim=mean(imstack,3);
    mim=filter2(h,mim);
    maxima=maximumfindcall(mim);
    figure(88);
    imagesc(mim);
    int=maxima(:,3);
    
    mimc=mim(p.roisize:end-p.roisize,p.roisize:end-p.roisize);
    mmed=myquantile(mimc(:),0.1);
    imt=mim(mim<mmed);
    cutoff=mean(imt(:))+3*std(imt(:));
    maxima=maxima(int>cutoff,:);
    hold on
    plot(maxima(:,1),maxima(:,2),'mo')
    hold off
    numframes=size(imstack,3);
    bind=length(b)+size(maxima,1);
    for l=1:size(maxima,1)
        b(bind).loc.frames=(1:numframes)';
        b(bind).loc.filenumber=zeros(numframes,1)+k;
        b(bind).filenumber=k;
        b(bind).pos=maxima(l,1:2);
        b(bind).stack.image=imstack(b(bind).pos(2)+rsr,b(bind).pos(1)+rsr,:);
        b(bind).stack.framerange=1:numframes;
        bind=bind-1;
    end
    fmax=max(fmax,numframes);
end

p.fminmax=[1 fmax];
p.cam_pixelsize_um=[1 1]/1000;
p.pathhere=fileparts(filelist{1});
end