function out=simulate2c
% transformfile='/Users/jonas/Documents/Data/ROI2_20per639_50msexp_1/ROI2_20per639_50msexp_1/ROI2_20per639_50msexp_1_MMStack2_T.mat';
cal3dfile='/Users/jonas/Documents/Data/ROI2_20per639_50msexp_1/ROI2_20per639_50msexp_1_3dcal.mat'; 
% transformfile='/Users/jonas/Documents/Data/ROI2_20per639_50msexp_1/simulstack_affine_2_T.mat';
transformfile='/Users/jonas/Documents/Data/ROI2_20per639_50msexp_1/affine_T.mat';
% transformfile='/Users/jonas/Documents/Data/ROI2_20per639_50msexp_1/testmirrorA_T.mat';
% cal3dfile='/Users/ries/Documents/Data/3D/global3D/ROI2_20per639_50msexp_1_3dcal.mat'; 
% transformfile='/Users/ries/Documents/Data/3D/global3D/ROI2_20per639_50msexp_1/ROI2_20per639_50msexp_1_MMStack_T.mat';

mirror=true;
 midp=255;
nbeads=25;
roisize=31;
dr=(roisize-1)/2;
zpos=1:201; %in slices
% zpos=99:101;
pixelsize=100;
Nphot=16000;
Nbg=5;
img=zeros(512,512,length(zpos),'single');

x=rand(nbeads,1)*(size(img,1)-roisize*2)+roisize;
y=rand(nbeads,1)*(size(img,2)-roisize*2)/2+roisize;

xnm=x*pixelsize; ynm=y*pixelsize;

t=load(transformfile); transformation=t.transformation;
t=load(cal3dfile); coeffh=t.SXY.cspline.coeff;


% coeff=coeffh;
coeff{1}=coeffh;coeff{2}=coeffh;
[xtnm,ytnm]=transformation.transformCoordinatesFwd(xnm,ynm);
xt=xtnm/pixelsize;yt=ytnm/pixelsize;
figure(88);plot(x,y,'+',xt,yt,'x')
% coord=horzcat(y,x);
coord=horzcat(x,y);

channel=1;

makeimg;
% coord=horzcat(yt,xt);
coord=horzcat(xt,yt);
channel=2;
makeimg
img=img+Nbg;


if mirror
   
    img(midp:end,:,:)=img(end:-1:midp,:,:);
    
end
imageslicer(img);

imgnoise=poissrnd(img);
outfile=[fileparts(cal3dfile) filesep 'simulstackM.tif'];
outfile2=[fileparts(cal3dfile) filesep 'simulstackMp.tif'];
% imgnoise=permute(imgnoise,[2,1,3]);
saveastiff(uint16(imgnoise),outfile);
saveastiff(uint16(img+Nbg),outfile2);
function makeimg
roipos=round(coord);
incoord=coord-roipos+dr;
incoord(:,[1 2])=incoord(:,[2 1]);
incoord(:,4)=Nphot;
incoord(:,5)=0;
for zk=1:length(zpos)
    incoord(:,3)=zpos(zk);
    impsf=renderPSF(coeff{channel},incoord,roisize);
    for k=1:nbeads
        if all(roipos(k,:)>dr)&&all(roipos(k,:)<size(img,1)-dr+1)
       img(roipos(k,2)-dr:roipos(k,2)+dr,roipos(k,1)-dr:roipos(k,1)+dr,zk)=img(roipos(k,2)-dr:roipos(k,2)+dr,roipos(k,1)-dr:roipos(k,1)+dr,zk)+impsf(:,:,k);
    
        end
    end
   
end
end
end