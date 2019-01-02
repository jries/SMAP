function radialfourier

file='/Users/jonas/OneDrive/Projects/NPC_Reference/images/Figure1/ExM/NUP96_WF.tif';
file='/Users/jonas/OneDrive/Projects/NPC_Reference/images/Figure1/sted.png';
img=imread(file);
img=double(img);
if ndims(img)>2
    img=sum(img,3);
end
imf=imfinfo(file);
pixelsize=1000/imf.XResolution %pixelsize in nm

in1=fftshift(fft2(img));
rs=radialsum(abs(in1));
figure(88)
subplot(2,2,1)
imagesc(img)
subplot(2,2,2);
imagesc((abs(in1)));
subplot(2,2,3)
rpix=(1:length(rs))';
qmax=0.5/pixelsize;

freq=linspace(0,qmax, length(rs));


plot(freq,(rs./rpix))
ax=gca;
xt=ax.XTick;
for k=1:length(xt)
    xtl{k}=num2str(1/xt(k),'%2.0f');
end
ax.XTickLabel=xtl;
ax.XTickMode='manual';
end
function rs=radialsum(img)
s=size(img);
center=floor((s+1)/2);
rs=zeros(ceil(s(1)/2)+1,1);
for k=1:s(1)
    for l=1:s(2)
        d=sqrt((k-center(1)).^2+(l-center(2)).^2);
        ind=round(d)+1;
        if ind<=length(rs)
        rs(ind)=rs(ind)+img(k,l);
        end
    end
end
end

