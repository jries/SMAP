function makeBlinkMovie(loc,filestruc,p,obj)
global reconstructgui_abort
gaussfac=0.7;
pix_cam=filestruc.info.cam_pixelsize_um*1000;
roi=filestruc.info.roi;
path=fileparts(filestruc.name); %rather in top class, pass on
[file, path]=uigetfile([path filesep '*.tif'],'Select single images');
df=myfastdir(path,  'img*.tif');
pfile.filelist=df;
pfile.path=path;
pfile.bufferSize=200;

hc=figure;


% loc.frame=loc.frame+1; %filtered wth numberingroup>1: take the next (brighrter?) for grouped data.

% tifloader=plugintemp.plugins('WorkflowModules','Loaders','TifLoader',hc,par);
% 
% getlocs=plugintemp.plugins('WorkflowModules','IntensityCaluclator','pushLocs',hc,par);
% getlocs.setInputModule(1,tifloader); 
% roicutter=makeplugin(obj,{'WorkflowModules','Peakfinders','RoiCutterWF'},hc,par);
% roicutter.setInputModule(1,tifloader); 
% roicutter.setInputModule(2,getlocs); 
fileloader=imageloaderAll([path file],[],obj.P);
% fileloader=fitter.FileLoader;
% fileloader.setParameters(pfile);

sr=round(2*pix_cam(1)/p.sr_pixrec);
% pcut.kernelSize=13;
% roiCutter=fitter.RoiCutter;
% roiCutter.setParameters(pcut);

rxpix1=(((p.sr_pos(1)-p.sr_size(1))/pix_cam(1):(p.sr_pos(1)+p.sr_size(1))/pix_cam(1))-roi(1));
rypix1=(((p.sr_pos(2)-p.sr_size(2))/pix_cam(2):(p.sr_pos(2)+p.sr_size(2))/pix_cam(2))-roi(2));
rxpix=round(rxpix1);
rypix=round(rypix1);

ddx=rxpix(1)-rxpix1(1);
ddy=rypix(1)-rypix1(1);
% rx=((p.sr_pos(1)-p.sr_size(1))/p.sr_pixrec:(p.sr_pos(1)+p.sr_size(1))/p.sr_pixrec);
% ry=((p.sr_pos(2)-p.sr_size(2))/p.sr_pixrec:(p.sr_pos(2)+p.sr_size(2))/p.sr_pixrec);
% 
% 
% rx=round(rx-roi(1)*pix_cam(1)/p.sr_pixrec);
% ry=round(ry-roi(2)*pix_cam(2)/p.sr_pixrec);

pos.x=loc.xnm;pos.y=loc.ynm;pos.s=max(loc.locprecnm*gaussfac,p.sr_pixrec/2);ch=loc.channel;
rangex=[p.sr_pos(1)-p.sr_size(1) p.sr_pos(1)+p.sr_size(1)];
rangey=[p.sr_pos(2)-p.sr_size(2) p.sr_pos(2)+p.sr_size(2)];

inch1=ch==1;
inch2=ch==2;
pos1.x=pos.x(inch1);pos1.y=pos.y(inch1);pos1.s=pos.s(inch1);
pos2.x=pos.x(inch2);pos2.y=pos.y(inch2);pos2.s=pos.s(inch2);

[srimfinal1,nlocs,G]=gaussrender(pos1,rangex, rangey, p.sr_pixrec, p.sr_pixrec);
[srimfinal2,nlocs,G]=gaussrender(pos2,rangex, rangey, p.sr_pixrec, p.sr_pixrec);

srimadd1=0*srimfinal1; %
srimadd2=0*srimfinal1; %
average=0;

ssr=size(srimfinal1);
ax1=initaxis(p.resultstabgroup,'running movie');
% ax2=initaxis(p.resultstabgroup,'final');
axes(ax1)

frames=p.frame_min+p.numberOfFrames:-1:p.frame_min;

frand=rand(length(frames),1);
[~,finds]=sort(frand);
frames=frames(finds);

% mov(length(frames))=struct('cdata',[],'colormap',[]);
% indmov=1;
ind=strfind(path,filesep);
p2=path(1:ind(end-1));
switch p.outputFormat.selection
    case {'MPEG-4'}
        ext='.mp4';
    case {'Uncompressed AVI','Motion JPEG 2000'}
        ext='.avi';
end
[f,path2]=uiputfile(['*' ext],'select movie', [p2 ext]);

outfile=[path2 f];
% aviobj=VideoWriter(outfile,'Archival');
aviobj=VideoWriter(outfile,p.outputFormat.selection);
% tifobj=Tiff(outfile,'w');
% 
% tagstruct.ImageLength = ceil(ssr(1)/16)*16;
% tagstruct.ImageWidth = ceil(ssr(2)*2/16)*16;
% tagstruct.Photometric = Tiff.Photometric.RGB;
% tagstruct.BitsPerSample = 8;
% tagstruct.SamplesPerPixel = 3;
% tagstruct.RowsPerStrip = 16;
% tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
% tagstruct.Software = 'MATLAB';
% tagstruct.Compression=7;
% tifobj.setTag(tagstruct)



% aviobj.Quality=100;
aviobj.FrameRate=p.framerate;
open(aviobj);

normim=[];
% tifstack=zeros(ssr(1),ssr(2)*2,3,length(frames));
for k=1:length(frames)
ind=find(loc.frame==frames(k));
if ~isempty(ind)
  
pos.x=loc.xnm(ind);pos.y=loc.ynm(ind);pos.s=max(loc.locprecnm*gaussfac,p.sr_pixrec);ch=loc.channel(ind);
inch1=ch==1 | ch==0;
inch2=ch==2;
pos1.x=pos.x(inch1);pos1.y=pos.y(inch1);pos1.s=pos.s(inch1);
pos2.x=pos.x(inch2);pos2.y=pos.y(inch2);pos2.s=pos.s(inch2);

srimhere1=gaussrender(pos1,rangex, rangey, p.sr_pixrec, p.sr_pixrec,[],[],G);
srimadd1=srimadd1+srimhere1;

srimhere2=gaussrender(pos2,rangex, rangey, p.sr_pixrec, p.sr_pixrec,[],[],G);
srimadd2=srimadd2+srimhere2;

imf=double(fileloader.getimage(frames(k)));
emccd=1;
if emccd
imf=imf(:,end:-1:1);
end
% imf=double(fileloader.getImage(k));
if ~isempty(imf)

    imcut=imf(rypix,rxpix);
%     xpix=(pos.x/pix_cam(1)-roi(1));
%     ypix=(pos.y/pix_cam(2)-roi(2));
    

imcut=imresize(imcut,pix_cam(1)/p.sr_pixrec,'nearest');
% imcut=imb(rx,ry);
% 
imcut=imcut-min(imcut(:));
if isempty(normim)
    normim=myquantile(imcut(:),0.99997);
end
normim2=myquantile(imcut(:),0.99997);
normim=max(0.8*normim+0.2*normim2,normim2);
imcut=imcut/normim;
% 
average=imcut+average;
% 
% 
xp=(loc.xnm(ind)-p.sr_pos(1)+p.sr_size(1))/p.sr_pixrec-((ddx-.5)*pix_cam(1)/p.sr_pixrec);
yp=(loc.ynm(ind)-p.sr_pos(2)+p.sr_size(2))/p.sr_pixrec-((ddy-.5)*pix_cam(2)/p.sr_pixrec);

ims=plotsquares(imcut,yp,xp,sr,ssr);
% 
% 
% 
% 
% 
% 
% 
collage=makecollage(srimadd1,ims,srimadd2);

collage=uint8(collage*255);
% 
writeVideo(aviobj,collage);
% % write(tifobj,uint8(collage*255));
% % indmov=indmov+1;
% 
image(collage,'Parent',ax1);
axis equal
% colormap hot
drawnow
% 
% 
%     if reconstructgui_abort
%         error('stopped manually')
%         
%     end
end
end
end

averageout=plotsquares(average,[],[],sr,ssr);
averageout=averageout-min(averageout(:));

% axes(ax2);
collage=makecollage(srimfinal1/max(srimfinal1(:))*3,averageout/myquantile(averageout(:),0.99997),srimfinal2/max(srimfinal2(:))*3);
% image(collage);
% figure(29);movie(mov)
% axis equal
% writeVideo(aviobj,im2frame(collage));
close(aviobj);
% tifobj.close();
imwrite(collage,[outfile(1:end-4) '_all.tif'])


function collage=makecollage(srimage,dlimage,srim2)
ssr=size(srimage);
norm=max(myquantile(srimage(:),0.9996));
if norm==0
    norm=max(1,max(srimage(:)));
end
srimplot=srimage/norm;
srimplot(srimplot>1)=1;
srimplot(1,1)=1;

norm2=max(myquantile(srim2(:),0.9994));
if norm2==0
    norm2=max(1,max(srim2(:)));
end
srimplot2=srim2/norm2;
srimplot2(srimplot2>1)=1;
srimplot2(1,1)=1;

lut=hot(257);
srimplotrgb=ind2rgb(ceil(srimplot*255+1),lut);
srimplotrgb2=0*srimplotrgb;
srimplotrgb2(:,:,2)=(srimplot2);
srimplotrgb2(:,:,3)=(srimplot2);

srimplotrgb=srimplotrgb+srimplotrgb2;
collage=zeros(ssr(1),ssr(2)*2,3);
for c=1:3
    collage(:,1:ssr(2),c)=dlimage(:,:,c);
    collage(:,ssr(2)+1:2*ssr(2),c)=srimplotrgb(:,:,c);
end
collage(collage>1)=1;

function imoutf=plotsquares(image,x,y,sr,ssr)
s=size(image);
badind=x<1|x>s(1)|y<1|y>s(2);
x(badind)=[];
y(badind)=[];
imout=zeros(s(1)+2*sr+4,s(2)+2*sr+4,3);
sim=size(image);
for k=1:3
    imout(sr+1:sr+sim(1),sr+1:sr+sim(2),k)=image;
end
x=round(x+sr);
y=round(y+sr);
sc=5;
for loc=1:length(x)
    imout(x(loc)-sr:x(loc)+sr,y(loc)+sr,1:3)=1;
    imout(x(loc)-sr:x(loc)+sr,y(loc)-sr,1:3)=1;
    imout(x(loc)+sr,y(loc)-sr:y(loc)+sr,1:3)=1;
    imout(x(loc)-sr,y(loc)-sr:y(loc)+sr,1:3)=1;
%     
%     imout(x(loc)-sc:x(loc)+sc,y(loc),:)=0;
%     imout(x(loc),y(loc)-sc:y(loc)+sc,:)=0;
%         imout(x(loc)-sc:x(loc)+sc,y(loc),1)=1;
%     imout(x(loc),y(loc)-sc:y(loc)+sc,1)=1;
end

imoutf=imout(sr+1:sr+ssr(1),sr+1:sr+ssr(2),:);