function [out,dxo,dyo]=myhist2(xg,yg,dx,dy,mx,my,w)
%x
%mx(1),mx(2),my minimum, maximum values
%dx pixelsize for reconstruction 

if nargin<7
    w=1+0*xg;
end
if nargin<6
    mx=[min(xg) max(xg)];
    my=[min(yg) max(yg)];
end

minx=floor(mx(1)/dx);
miny=floor(my(1)/dy);
maxx=ceil(mx(2)/dx);
maxy=ceil(my(2)/dy);

x=round(xg/dx);
y=round(yg/dy);



x=x-(minx)+1;y=y-(miny)+1;

maxx2=((maxx)-minx);maxy2=(maxy-miny);
maxim=1e9;

if maxx2*maxy2<maxim
    sw=size(w);
  out=zeros(maxx2,maxy2,sw(2));
else
    disp([ 'reconstruction image too big > ' num2str(maxim)])
    asf
end
sr=size(out);
ig=( x>0&y>0&x<sr(1)&y<sr(2));    
fig=find(ig);
so=size(out);

for c=1:sw(2)
        linind=sub2ind(so,x(ig),y(ig),c+0*fig);
        for k=1:length(linind)
            out(linind(k))=out(linind(k))+w(fig(k),c);
        end
end
% out(1,1)=0;
[m,mi]=max(out(:)); out(mi)=out(mi+1);
% imagesc([mx(1) mx(2)],[my(1) my(2)],out');

dxo=(minx:maxx-1)*dx;
dyo=(miny:maxy-1)*dy;
    
out(:,1,:)=[];
out(1,:,:)=[];


