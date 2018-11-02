function imout=gettif(file)
imout.image=imread(file); 
sim=size(imout.image);
imout.info.Width=sim(1);
imout.info.Height=sim(2);
imout.info.roi=getRoiTif(file);
if all((imout.info.roi(:))==([0 ;0 ;512; 512]))
    imout.info.roi(3:4)=sim;
end
imout.info.name=file;
end