function [out,dxo,dyo,dzo]=myhist3(xg,yg,zg,d,mx,my,mz,w)
%x
%mx(1),mx(2),my minimum, maximum values
%dx pixelsize for reconstruction 

if nargin<8
    w=ones(length(xg),1);
end
if nargin<7
    mx=[min(xg) max(xg)];
    my=[min(yg) max(yg)];
    mz=[min(zg) max(zg)];
end
if length(d)<3
    d(2)=d(1);d(3)=d(1);
end
minx=floor(mx(1)/d(1));
miny=floor(my(1)/d(2));
minz=floor(mz(1)/d(3));
maxx=ceil(mx(2)/d(1));
maxy=ceil(my(2)/d(2));
maxz=ceil(mz(2)/d(3));

x=round(xg/d(1));
y=round(yg/d(2));
z=round(zg/d(3));


x=x-(minx)+1;y=y-(miny)+1;z=z-(minz)+1;

maxx2=((maxx)-minx);maxy2=(maxy-miny);maxz2=(maxz-minz);
maxim=1e9;

if maxx2*maxy2*maxz2<maxim
    sw=size(w);
  out=zeros(maxx2,maxy2,maxz2,sw(2));
else
    disp([ 'reconstruction image too big > ' num2str(maxim)])
    asf
end
sr=size(out);
ig=( x>0&y>0&z>0&x<=sr(1)&y<=sr(2)&z<=sr(3));    
fig=find(ig);
so=size(out);

for c=1:sw(2)
        linind=sub2ind(so,x(ig),y(ig),z(ig),c+0*fig);
        for k=1:length(linind)
            out(linind(k))=out(linind(k))+w(fig(k),c);
        end
end
% out(1,1)=0;
[m,mi]=max(out(:)); out(mi)=out(mi+1); %?????
% imagesc([mx(1) mx(2)],[my(1) my(2)],out');

dxo=(minx:maxx-1)*d(1);
dyo=(miny:maxy-1)*d(2);
  dzo=(minz:maxz-1)*d(3);
  
out(:,1,:)=[];
out(1,:,:)=[];
out(:,:,1)=[];

