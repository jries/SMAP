%overlay dual color tracking images

%obj=g;
datafile='/Users/ries/Data_local/temp/kinesins/614_purified_DualChannel_488i3_640i5_Alternate_3/614_purified_DualChannel_488i3_640i5_Alternate_MMStack_Pos0.ome.tif';
Tfile='/Users/ries/Data_local/temp/kinesins/Bead_Cal_488_640_lessIntensityLaser_T.mat';

alternating=true;

tt=load(Tfile).transformation;
il=imageloaderMM; il.attachPar(g.P);
il.openi(datafile);
il.prefit;

numf=il.metadata.numberOfFrames;

sx=ceil(il.metadata.Width/2);
imcomb=zeros(il.metadata.Height,sx,3,numf);


for k=1:numf
    img=double(il.getimage(k));
    imgt=tt.transformImageToTarget(2,img,'pixel',il.metadata.roi);
    imcomb(:,:,2,k)=img(:,1:sx);
    % imcomb(:,:,3,k)=img(:,1:sx);
    imcomb(:,:,1,k)=imgt(:,1:sx);
    
end

ih=imx(imcomb);

if alternating
    inodd=1:2:numf;
    ineven=2:2:numf;
    im1=imcomb(:,:,1,inodd);
    im2=imcomb(:,:,2,ineven);
    im3=imcomb(:,:,3,inodd);
    imcomb=zeros(size(imcomb,1),size(imcomb,2),3,ceil(numf/2));
    imcomb(:,:,1,:)=im1;imcomb(:,:,2,:)=im2;imcomb(:,:,3,:)=im3;
end

fout=strrep(datafile,'.ome.tif','_combined.tif')
options.color=true;
saveastiff(squeeze(single(imcomb)),fout,options);
% fout2=strrep(datafile,'.ome.tif','_combined2.tif');
% saveastiff(squeeze(single(imcomb(:,:,2,:))),fout2);