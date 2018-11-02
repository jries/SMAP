%mammalian cells: temporal reconstruction
%parameters
p.cutoff=4; % for isosurface, at max(V(:))/cutoff
p.cutoffdc=3; % for isosurface, at max(V(:))/cutoff
p.sgauss=2;%smoothing of volume
p.sgauss2=.6;%smoothing of volume
p.winsize=200; %plotting
p.minx=0; %x: use 0 for central cut
p.pixelsize=5; %nm, for volume reconstruction
% p.limxy=[0 40]; % for plotting, in volume pixels
% p.limz=[0 40];

p.polarfactor=[.5 1]; %intensity factor for layer1, layer2
rot=1; % rotate sites to align bottom holes

timewindows=20; % number of sites per average: (number of sites / this number)
timepoints=100; % increment: (number of sites / this number)
framerate=25;

global se
sites=[];
use=[];
use=getFieldAsVector(se.sites,'annotation','use');
sites=se.sites(use);


if ~isfield(sites(1).evaluation.CME3DDSpherefit.map3D,'coordinates2')
    for k=1:length(sites)
        sites(k).evaluation.CME3DDSpherefit.map3D.coordinates2=struct('x',[],'xc',[],'y',[],'yc',[],'z',[],'zc',[]);
    end
end

% rSphere=getFieldAsVector(sites,'evaluation','CME3DDSpherefit','map3D','rSphere');
mainFraction=getFieldAsVector(sites,'evaluation','CME3DDSpherefit','map3D','mainFraction');

[~,indsort]=sort(mainFraction);

sitessort=sites(indsort);

% sitessort=sites;
% sitessort(1:10)=[];
rSphere=getFieldAsVector(sitessort,'evaluation','CME3DDSpherefit','map3D','rSphere');
mainFraction=getFieldAsVector(sitessort,'evaluation','CME3DDSpherefit','map3D','mainFraction');
topcoverage=getFieldAsVector(sitessort,'evaluation','CME3DDSpherefit','map3D','topcoverage');
thetac=getFieldAsVector(sitessort,'evaluation','CME3DDSpherefit','map3D','thetac');
phic=getFieldAsVector(sitessort,'evaluation','CME3DDSpherefit','map3D','phic');
thetacoverage=getFieldAsVector(sitessort,'evaluation','CME3DDSpherefit','map3D','bwfit','thetacoverage');

dN=round(length(indsort)/timewindows);
df=round(length(indsort)/timepoints);

range=1:df:length(indsort)-dN+df;
range(end)=min(range(end),length(sitessort)-dN);

rmean=zeros(1,length(indsort));
rcorr=zeros(1,length(indsort));
F(length(range))=struct('cdata',[],'colormap',[]);
polall=[];
for k=1:length(range)
    numloc=0;
    numloc2=0;
    
    indh=range(k):range(k)+dN;
    
    for s=indh
        numloc=numloc+length(sitessort(s).evaluation.CME3DDSpherefit.map3D.coordinates.r);  
        numloc2=numloc2+length(sitessort(s).evaluation.CME3DDSpherefit.map3D.coordinates2.x); 
    end
    
        
    thcov=mean(pi-thetacoverage(indh));
    Rhere=mean(rSphere(indh));
        rmean(indh)=Rhere;

    rcorr(indh)=rmean(indh)./rSphere(indh);   

    x=zeros(1,numloc);y=zeros(1,numloc);z=zeros(1,numloc);
    x2=zeros(1,numloc2);y2=zeros(1,numloc2);z2=zeros(1,numloc2);
    
    ind=1;ind2=1;
    for s=indh
%         if topcoverage(s)
%             tcorr=0;
%         else
%             tcorr=pi;
%         end
%         pcorr=0;
% %         tcorr=-thetac(s);pcorr=-phic(s);
        if rot
                xh=sitessort(s).evaluation.CME3DDSpherefit.map3D.coordinates.xc;
                x(ind:ind+length(xh)-1)=xh*rcorr(s);
                y(ind:ind+length(xh)-1)=sitessort(s).evaluation.CME3DDSpherefit.map3D.coordinates.yc*rcorr(s);
                z(ind:ind+length(xh)-1)=sitessort(s).evaluation.CME3DDSpherefit.map3D.coordinates.zc*rcorr(s);

                    xh2=sitessort(s).evaluation.CME3DDSpherefit.map3D.coordinates2.xc;
                    x2(ind2:ind2+length(xh2)-1)=xh2*rcorr(s);
                    y2(ind2:ind2+length(xh2)-1)=sitessort(s).evaluation.CME3DDSpherefit.map3D.coordinates2.yc*rcorr(s);
                    z2(ind2:ind2+length(xh2)-1)=sitessort(s).evaluation.CME3DDSpherefit.map3D.coordinates2.zc*rcorr(s);
        %         x(ind:ind+length(r)-1)=r;
        %         rch(ind:ind+length(r)-1)=r*rcorr(s);
        %         th(ind:ind+length(r)-1)=sitessort(s).evaluation.CME3DDSpherefit.map3D.coordinates.t+tcorr;
        %         ph(ind:ind+length(r)-1)=sitessort(s).evaluation.CME3DDSpherefit.map3D.coordinates.p+pcorr;

        else
                xh=sitessort(s).evaluation.CME3DDSpherefit.map3D.coordinates.x;
                x(ind:ind+length(xh)-1)=xh*rcorr(s);
                y(ind:ind+length(xh)-1)=sitessort(s).evaluation.CME3DDSpherefit.map3D.coordinates.y*rcorr(s);
                z(ind:ind+length(xh)-1)=sitessort(s).evaluation.CME3DDSpherefit.map3D.coordinates.z*rcorr(s);
                    xh2=sitessort(s).evaluation.CME3DDSpherefit.map3D.coordinates2.x;
                    x2(ind2:ind2+length(xh2)-1)=xh2*rcorr(s);
                    y2(ind2:ind2+length(xh2)-1)=sitessort(s).evaluation.CME3DDSpherefit.map3D.coordinates2.y*rcorr(s);
                    z2(ind2:ind2+length(xh2)-1)=sitessort(s).evaluation.CME3DDSpherefit.map3D.coordinates2.z*rcorr(s);        
        end        
        ind=ind+length(xh);
        ind2=ind2+length(xh2);
        
    end
    [f,imcomp]=plotcoords(x,y,z,Rhere,thcov,p,x2,y2,z2);
    if isempty(polall)
        s=size(imcomp);
        polall=zeros(s(1),s(2),s(3),length(range));
    end
    polall(:,:,:,k)=imcomp;
   
    F(k)=getframe(f);
end

temp=polall(:,:,1,:);
m1=quantile(temp(:),0.998);
temp=polall(:,:,2,:);
m2=quantile(temp(:),0.999);
polalln=polall;
polalln(:,:,1,:)=polalln(:,:,1,:)/m1;
polalln(:,:,2,:)=polalln(:,:,2,:)/m2;
%polalln(polalln>1)=1;

h=fspecial('gauss',5,1.2);
if 1 %if you want to save polar plot
    for k=1:size(polall,4)
        imh=polalln(:,:,:,k);
        imhf=imfilter(imh,h);
        f=figure(111);
        imagesc([-1 1]*p.winsize,[-1 1]*p.winsize,imhf);
        drawnow
        
        F(k)=getframe(f);
    end
end

path=fileparts(se.files(1).name);
[f,path]=uiputfile([path filesep 'movie.mp4']);
if f
    v=VideoWriter([path f],'MPEG-4');
    v.FrameRate=framerate;
    open(v)
    for k=1:length(F)
        writeVideo(v,F(k));
    end
    close(v);
end
%f=figure(123);movie(f,F,2)



function [f,imcomp]=plotcoords(x,y,z,R,thcov,p,x2,y2,z2)
% figure(88)
% scatter3(x,y,z+R,[],z+R)
% axis equal
% [t,p,r]=cart2sph(x,y,z);
% figure(89)
% polarplot(p,r,'.')
dz=R*cos(thcov);
dz=R;
xrange=[-p.minx,p.winsize];
yrange=[-p.winsize,p.winsize];
zrange=[-p.winsize,p.winsize]-p.winsize/2;

xr=xrange(1)+p.pixelsize:p.pixelsize:xrange(2)-p.pixelsize;
yr=yrange(1)+p.pixelsize:p.pixelsize:yrange(2)-p.pixelsize;
zr=zrange(1)+p.pixelsize:p.pixelsize:zrange(2)-p.pixelsize;
V=myhist3(x,y,-(z+dz),p.pixelsize,xrange,yrange,zrange);
Vs=smooth3(V,'gaussian',[1 1 1]*(round(2*p.sgauss/2)*2+1),p.sgauss);


V2=myhist3(x2,y2,-(z2+dz),p.pixelsize,xrange,yrange,zrange);
Vs2=smooth3(V2,'gaussian',[1 1 1]*(round(2*p.sgauss/2)*2+1),p.sgauss2);

[X,Y,Z]=meshgrid(yr,xr,zr);

co=max(Vs(:))/p.cutoff;
codc=max(Vs(:))/p.cutoffdc;

f1=figure(90);
clf
hiso=patch(isosurface(X,Y,Z,Vs,co),'EdgeColor','none','FaceColor',[1,.75,.65]);
isonormals(X,Y,Z,Vs,hiso);

hcap=patch(isocaps(X,Y,Z,Vs,co),'FaceColor','interp','EdgeColor','none');
axis equal
view(-35,30)

xlim(yrange)
ylim(xrange)
zlim(zrange)

lightangle(45,30);
lighting gouraud

f2=figure(91);
clf
hiso=patch(isosurface(X,Y,Z,Vs,codc),'EdgeColor','none','FaceColor',[1,.75,.65],'FaceAlpha',0.75);
isonormals(X,Y,Z,Vs,hiso);
hold on
if 0
indin=x2>0;
scatter3(y2(indin),x2(indin),-z2(indin)-dz,2,'k')
else
codc2=max(Vs2(:))/p.cutoffdc;
% codc2=.5;
hiso2=patch(isosurface(X,Y,Z,Vs2,codc2),'EdgeColor','none','FaceColor',[0,.75,1]);
isonormals(X,Y,Z,Vs2,hiso2);
hcap2=patch(isocaps(X,Y,Z,Vs2,codc2),'FaceColor','interp','EdgeColor','none');
end

axis equal
view(-35,30)

xlim(yrange)
ylim(xrange)
zlim(zrange)

lightangle(45,30);
lighting gouraud



%2D analysis

img=getpolarimage(x,y,z,R,p)*p.polarfactor(1);
img2=getpolarimage(x2,y2,z2,R,p)*p.polarfactor(2);
s=size(img);
imcomp=zeros(s(2),s(1),3);
imcomp(:,:,1)=img';
imcomp(:,:,2)=img2';
f3=figure(45);
imagesc(imcomp);
axis equal

f=f1; %select which figure is saved
% indg=r~=0;
% % [z,x]=pol2cart(p(indg),r(indg));
% [z,x,y]=sph2cart(t(indg),p(indg),r(indg));
% zm=z+mean(r);
% figure(87);
% plot(x,-zm,'.')
% % xlim([-200 200])
% % ylim([-200 200])
% % posp.x=-zm;posp.y=x;
% ps=5;
% edgesx=-200:ps:200;
% edgesz=-100:ps:300;
% srim=histcounts2(x,zm,edgesx,edgesz);
% [Y,X]=meshgrid(rangez(1)+pixelsz/2:pixelsz:rangez(2),ranger(1)+pixelsx/2:pixelsx:ranger(2));

% srim=histrender(posp,[-200 200], [-200 200], 5, 5);
% figure(89)
% imagesc(srim')
% axis('equal')

% % subplot(2,3,1)
% 
% imn=srim./(X'+pixelsx/4);
% % imn=srim';
% imn2=[imn(:,end) imn(:,end:-1:2) imn];

% figure(88);polarplot(th,rch,'.')
% rlim([0,max(rmean)*1.1])
drawnow;
% waitforbuttonpress
% pause(.2)

end



function imgo=getpolarimage(x,y,z,R,p)
[a,e,r]=cart2sph(x,y,z);
[xp,zp]=pol2cart(e,r);
zp2=zp+R;
w=1./(abs(cos(e))+0.2)/R;
img=myhist2(xp,zp2,p.pixelsize,p.pixelsize,[0 1]*p.winsize,[-1 1]*p.winsize,w');
imgo=vertcat(img(end:-1:2,:),img);
end