%overlay dual color tracking images

obj=g;
datafile='/Users/ries/datalocal/2color_kinesin/25_50ms_561nm01_640nm02_600w52_676w37_1_MMStack_Default.ome.tif';
Tfile='/Users/ries/datalocal/2color_kinesin/75_50ms_561nm01_640nm02_600w52_676w37_1_T.mat';

tt=load(Tfile).transformation;
il=imageloaderMM; il.attachPar(g.P);
il.openi(datafile);
il.prefit;

numf=il.metadata.numberOfFrames+1;

sx=ceil(il.metadata.Width/2);
imcomb=zeros(il.metadata.Height,sx,3,numf);


for k=1:numf
    img=double(il.getimage(k));
    imgt=tt.transformImageToTarget(2,img,'pixel',il.metadata.roi);
    imcomb(:,:,1,k)=img(:,1:sx);
    % imcomb(:,:,3,k)=img(:,1:sx);
    imcomb(:,:,2,k)=imgt(:,1:sx);
    
end

ih=imx(imcomb);

fout=strrep(datafile,'.ome.tif','_combined.tif')
options.color=true;
saveastiff(squeeze(single(imcomb)),fout,options);
% fout2=strrep(datafile,'.ome.tif','_combined2.tif');
% saveastiff(squeeze(single(imcomb(:,:,2,:))),fout2);