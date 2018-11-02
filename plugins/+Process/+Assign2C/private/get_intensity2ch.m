function loco=get_intensity2ch(loc,p,ind)

%parameters
if nargin<3
    ind=true(length(loc.(p.assignfield1.selection)),1);
end
if isnumeric(p.split_slope)
    slope=p.split_slope;
else
    slope=str2num(p.split_slope);
end
if length(slope)==1;
    slope(2)=slope(1);
end


if isnumeric(p.split_offset)
    offset=p.split_offset;
else
    offset=str2double(p.split_offset);
end

if isnumeric(p.split_edge)
    se=p.split_edge;
else
se=str2num(p.split_edge);
end

edge1=se(1);
if length(se)>1
    edge2=se(2);
else
    edge2=edge1;
end

if isnumeric(p.split_intmin)
    imin=p.split_intmin;
else
imin=str2num(p.split_intmin);
end
int1min=imin(1);
if length(imin)>1
    int2min=imin(2);
else
    int2min=int1min;
end

loco.channel=0*loc.channel+5;
int1=loc.(p.assignfield1.selection);
int2=loc.(p.assignfield2.selection);
% int1=loc.intA1;
% int2=loc.intB1;

if p.logscale
    int1(int1<1)=1;int2(int2<1)=1;
    int1=log10(int1);
    int2=log10(int2);
end
m1=myquantilefast(int1(int1~=0),[0.01,0.995],1e5);m2=myquantilefast(int2(int2~=0),[0.01,0.995],1e5);

map=max(m1(2), m2(2))+1;mip=min(m1(1),m2(1));

npix=250;

% ps=((map-mip)/npix);
ps1=(m1(2)-m1(1))/npix;
ps2=(m2(2)-m2(1))/npix;
% indgood=int1~=0&int2~=0;
% indgood=true(size(int1));
indgood=ind;
% img=myhist2(int1(indgood),int2(indgood),ps,ps,[mip map],[mip map]);
img=myhist2(int1(indgood),int2(indgood),ps1,ps2,m1,m2);
% img=myhist2(intf1,int2,ps,ps,[mip map],[mip map]);

ax1=initaxis(p.resultstabgroup,'i1 vs i2');
imagesc(m2,m1,img)
xlabel(p.assignfield2.selection)
ylabel(p.assignfield1.selection)
hold on
plotboundary
hold off

ax2=initaxis(p.resultstabgroup,'log i1 vs i2');
limg=log(img);
limg(isinf(limg))=-1;
imagesc(m2,m1,limg)
hold on
plotboundary
hold off
xlabel(p.assignfield2.selection)
ylabel(p.assignfield1.selection)
drawnow
i1co=slope(2)*int2+edge1+offset(1);
c1=int1>=i1co&int1>=int1min;
loco.channel(c1)=1;


i2co=1./slope(1)*((int1)+edge2-offset(1));
c2=int2>=i2co&int2>=int2min;
loco.channel(c2)=2;

loco.channel(int2==-100)=3;
loco.channel(int1==-100)=4;

ax3=initaxis(p.resultstabgroup,'log split');
imgc2=myhist2(int1(c2&ind),int2(c2&ind),ps1,ps2,m1,m2);
sout=size(imgc2);
outrgb=zeros(sout(1),sout(2),3);
outrgb(:,:,1)=imgc2;
imgc1=myhist2(int1(c1&ind),int2(c1&ind),ps1,ps2,m1,m2);
outrgb(:,:,2)=imgc1;

c3=~(c1|c2);
imgc3=myhist2(int1(c3&ind),int2(c3&ind),ps1,ps2,m1,m2);
outrgb(:,:,3)=imgc3;
outrgb=log(outrgb)+1;
outrgb(isinf(outrgb))=0;
imagesc(m2,m1,outrgb/myquantilefast(outrgb(:),.9995,1e5))
hold on
plotboundary
xlabel(p.assignfield2.selection)
ylabel(p.assignfield1.selection)
hold off
drawnow

function plotboundary
plot([mip map],slope(1)*[mip map]+offset(1),'w')
plot([mip map],slope(2)*[mip map]+offset(1),'w')
plot([mip map],slope(2)*[mip map]+edge1+offset(1),'m')
plot([mip map],slope(1)*[mip map]-edge2+offset(1),'m')
plot([mip 1/slope(2)*(int1min-edge1-offset(1))],(int1min)*[1 1],'c')
plot(int2min*[1 1],[mip slope(1)*int2min-edge2+offset(1)],'c')
end
end