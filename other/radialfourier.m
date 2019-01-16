function radialfourier
global path
[file, path]=uigetfile([path  '*.tif'],'MultiSelect','on');
if ~iscell(file)
    file={file};
end
f=figure(189)
hold off
for k=1:length(file)
    plotfunction([path file{k}])
end
legend(file)
end

function plotfunction(file)
img=imread([file]);
img=double(img);
if ndims(img)>2
    img=sum(img,3);
end
imf=imfinfo(file);
pixelsize=1000/imf.XResolution %pixelsize in nm

if contains(file,'MAX2_U2OS_Nup96GFP195_Wide_5i488_100nm_10slices_40ms_Roi2_2_1_MMStack_Pos0')
pixelsize= 131 %WF
% img=imresize(img,2);pixelsize=pixelsize/2;
end


% img=img-quantile(img,0.05); %no difference
in1=fftshift(fft2(img));
[rs,norm]=(radialsum(abs(in1).^2));
figure(188)
subplot(2,2,1)
imagesc(img)
subplot(2,2,2);
imagesc((abs(in1)));
f=figure(189);
f.Renderer='painters';

rpix=(1:length(rs))';
qmax=0.5/pixelsize;

freq=linspace(0,qmax, length(rs));

if contains(file,'STED')
  noise=mean(rs(end-25:5));  
else
  noise=mean(rs(end-2:end));
end

profnorm=(rs-rs(end))./norm;

% profnorm=profnorm-profnorm(end);
ind=find(freq>1/1500,1,'first');

profnorm=profnorm/profnorm(ind);
% plot(freq,(profnorm))
plot(freq,log(profnorm))
hold on;
ax=gca;
xt=ax.XTick;
xt=0:0.001:max(freq);
ax.XTick=xt;
for k=1:length(xt)
    xtl{k}=[num2str(1/xt(k),'% 2.0f')];
%     xtl{k}=['av';'cd'];
end
ax.XTickLabel=xtl;
ax.XTickMode='manual';
xlabel('1/frequency (nm)')
ylabel('power spectrum log(|F(I)|^2)')
 ylim([-10 1])
 xlim([0 0.025])    
end
function [rs,norm]=radialsum(img)
s=size(img);
center=floor((s+1)/2);
rs=zeros(ceil(s(1)/2)+1,1);
norm=0*rs;
for k=1:s(1)
    for l=1:s(2)
        d=sqrt((k-center(1)).^2+(l-center(2)).^2);
        ind=round(d)+1;
        if ind<=length(rs)
            rs(ind)=rs(ind)+img(k,l);
            norm(ind)=norm(ind)+1;
        end
    end
end
end

