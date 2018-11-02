function l=MakeNPCCoordinates
% radius=50;
% lpc=4; %number of corners
% dth=pi/48; %shift between corners
% theta=0:pi/4:2*pi-pi/1000; %8 corners
% thetaa=[];
% for k=1:lpc
%     thetaa=[thetaa theta+k*dth];
% end
% [l.x,l.y]=pol2cart(thetaa,radius);
% l.channel=ones(size(l.x));

corners=8;
%NUP 96
rN(1)=43.2;
rC(1)=47.2;
rN(2)=49.1;
rC(2)=52.5;

%total angle between corresponding points as in excel sheet, deg

aN(1)=9.1;
aC(1)=15.9;
aN(2)=15;
aC(2)=7.8;

r=rC;
a=aC;
%to deg, angle from ref:
a=a/2/180*pi;

r(3)=r(1);r(4)=r(2);
a(3)=-a(1);a(4)=-a(2);
z(1)=0;z(2)=50;z(3)=50;z(4)=0;

rall=zeros(length(r),corners);
angall=zeros(length(r),corners);
zall=zeros(length(r),corners);

dphi=0:2*pi/corners:2*pi-pi/1000;

for k=1:length(r)
    rall(k,:)=r(k);
    angall(k,:)=a(k)+dphi;
    zall(k,:)=z(k);
end
[l.x,l.y]=pol2cart(angall(:),rall(:));
l.z=zall(:);
l.channel=ones(size(l.x));