global se

savepos=0;
saveID=0;
saveinfo=0;
saveannot=0;
saveevaluation=0;
savename=0;
saveimage=0;
sum_image=se.sites(1).image.image;

% create summed and averaged image
for k=2:length(se.sites)
    sum_image=sum_image+se.sites(k).image.image;
end
av_image=sum_image/length(se.sites);
% figure(123);
% clf;
% imshow(rgb2gray(av_image(50:450,50:450,1:3)),[]);
% colormap(hot);



savepos={se.sites.pos};
saveID={se.sites.ID};
saveinfo={se.sites.info};
saveannot={se.sites.annotation};
saveevaluation={se.sites.evaluation};
savename={se.sites.name};
saveimage={se.sites.image};
numberofsites=length(se.sites);
% saveevaluation={se};
[FileName,PathName,FilterIndex]=uiputfile('*.mat','Save Evaluation Data As');
save(strcat(PathName,FileName),'savepos','saveID','saveinfo','saveannot','saveevaluation','savename','saveimage','av_image','sum_image','-v7.3');
save(strcat(PathName,FileName(1:end-4),'_images.mat'),'av_image','sum_image','numberofsites','-v7.3');

