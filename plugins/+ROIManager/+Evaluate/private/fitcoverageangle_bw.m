function fitp=fitcoverageangle_bw(bwim,startp,ax)
bwim=2*bwim-1;
sizein=size(bwim);
% dangle(1)=2*pi/sin(1);
% dangle(2)=pi/sin(2);
% % dangle=pi/64;
% tr=-pi/2:dangle(2):pi/2;
% pr=-pi:dangle(1):pi;
% r=linspace(-pi/2,pi/2,sin(2));
elsin=linspace(-1,1,sizein(1));
el=asin(elsin);
azimuth=linspace(-pi,pi,sizein(1));

[Elevation,Azimuth]=meshgrid(el,azimuth);

if nargin==1 || isempty(startp)
    startph1=[0,0,pi/2];
%     startph1=[pi/4,0,pi/2];
    imstart=coverage_sphere(startph1(1),startph1(2),startph1(3),Elevation,Azimuth);
    err1=sum((imstart(:)-bwim(:)).^2);
    startph2=[pi,0,pi/2];
    imstart=coverage_sphere(startph2(1),startph2(2),startph2(3),Elevation,Azimuth);
    err2=sum((imstart(:)-bwim(:)).^2);
    if err2<err1
        startph=startph2;
    else
        startph=startph1;
    end
else
    startph=startp;
    imstart=coverage_sphere(startph(1),startph(2),startph(3),Elevation,Azimuth);
end
options=optimset('lsqnonlin');
% options.Algorithm='levenberg-marquardt';


[fitp,resnorm]=lsqnonlin(@coverageerr,startph,[],[],options,Elevation,Azimuth,bwim);

% [fitp2,resnorm2]=lsqnonlin(@coverageerr,startph,[],[],options,Theta,Phi,-bwim);

% if resnorm2<resnorm
%     fitp=fitp2;
% end
%put back into range
th=mod(fitp(1),2*pi);
if th>pi
    th=th-2*pi;
end
ph=fitp(2);

if th>pi/2
    th=pi-th;
    ph=ph+pi;
elseif th<-pi/2
    th=-th-pi;
    ph=ph-pi;
end
ph=mod(ph,2*pi);
if ph>pi
    ph=ph-2*pi;
end
% [x,y,z]=sph2cart(fitp(2),fitp(1),1);
% [ph,th]=cart2sph(x,y,z);
 fitp(1)=th;fitp(2)=ph; 
    

if nargin>2
    hold(ax,'off')
imfit=coverage_sphere(fitp(1),fitp(2),fitp(3),Elevation,Azimuth);
% im=coverage(-pi/4,pi/4,pi/16,Theta,Phi);
imcombine=zeros(size(imstart,2),size(imstart,1),3);
imcombine(:,:,1)=bwim'*3;
imcombine(:,:,2)=imfit';
imcombine(:,:,3)=imstart';
imagesc(ax,azimuth,el,imcombine);
hold(ax,'on')
plot(ax,fitp(2)+pi/2,sin(fitp(1)),'*')
plot(ax,fitp(2)+pi/2-2*pi,sin(fitp(1)),'*')
end
% figure(81);imagesc(tr,pr,imfit);
end



function err=coverageerr(fitp,Theta,Phi,bwim)
% fitp
% fitp(3)=pi/4;
imtest=coverage_sphere(fitp(1),fitp(2),fitp(3),Theta,Phi);
err=imtest(:)-bwim(:);
if 0
imcombine=zeros(size(bwim,1),size(bwim,2),3);
imcombine(:,:,1)=(bwim+1)/2;
imcombine(:,:,2)=(imtest+1)/2;
% imcombine(:,:,3)=imstart;
figure(77);imagesc(imcombine);title(fitp)
drawnow
end
end