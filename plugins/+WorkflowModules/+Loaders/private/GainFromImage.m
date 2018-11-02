% function GainFromImage
img=double(img1);
dimg=double(img)-double(img2);
q=myquantilefast(img,[0.15,.2,.9]);
dI=q(2)-q(1);
int=min(img(:)):dI:q(3);
lint=length(int);
mi=zeros(lint-1,1);
vari=zeros(lint-1,1);
ds=250;
ix=1;
vart=0;
mit=0;
for k=1:lint-1
    ind=img>=int(k)&img<int(k+1);
    dimh=dimg(ind);
    imh=img(ind);
    indh=1:ds:length(imh);
    varih=inf;
    for l=1:length(indh)-1
%         varih=min(varih,));
        vart(ix)=var(dimh(indh(l):indh(l+1)));
        mit(ix)=mean(imh(indh(l):indh(l+1)));
        ix=ix+1;
    end
    vari(k)=varih;
%     vari(k)=var(dimg(ind))
    mi(k)=mean(img(ind));
end
indout=isinf(vari);
mi(indout)=[];vari(indout)=[];
figure(88)
plot(mit,vart/2,'.')
hold on
pix2adu=15.8/200;
plot(mit,(mit-offs)/pix2adu)
hold off

    
    
    