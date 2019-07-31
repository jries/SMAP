function fitp=fitcoverageangle_image(imin,startp,ax,rangex,rangey)
sin=size(imin);
% imin=sqrt(imin);
% dangle(1)=2*pi/sin(1);
% dangle(2)=pi/sin(2);
% % dangle=pi/64;
% tr=-pi/2:dangle(2):pi/2;
% pr=-pi:dangle(1):pi;
% r=linspace(-pi/2,pi/2,sin(2));
% trsin=linspace(-1,1,sin(1));
% tr=asin(trsin);
% pr=linspace(-pi,pi,sin(1));

trsin=double(linspace(rangey(1),rangey(2),sin(2)));
tr=asin(trsin);
pr=double(linspace(rangex(1),rangex(2),sin(1)));

[Theta,Phi]=meshgrid(tr,pr);

if length(startp)<4
    startp(4)=max(imin(:));%amplitude
%     startp(5)=0;%offset
end
if any(isnan(startp))
    fitp=[0 0 0];
    return
end
% if nargin==1 || isempty(startp)
%     startph1=[pi/4,0,pi/2];
%     imstart=coverage_sphere(startph1(1),startph1(2),startph1(3),Theta,Phi);
%     err1=sum((imstart(:)-bwim(:)).^2);
%     startph2=[pi/4+pi,0,pi/2];
%     imstart=coverage_sphere(startph2(1),startph2(2),startph2(3),Theta,Phi);
%     err2=sum((imstart(:)-bwim(:)).^2);
%     if err2<err1
%         startph=startph2;
%     else
%         startph=startph1;
%     end
% else
%     startph=startp;
% end
options=optimset('lsqnonlin');
% options.Algorithm='levenberg-marquardt';

imstart=startp(4)*(coverage_sphere(startp(1),startp(2),startp(3),Theta,Phi)+1);
try
[fitp,resnorm]=lsqnonlin(@coverageerr,startp,[],[],options,Theta,Phi,imin);
catch err
    startp
    err
    fitp=[0 0 0];
end

% fitp(1)=fitp(1)-rangey(1);
% fitp(2)=fitp(2)-rangex(1);
% [fitp2,resnorm2]=lsqnonlin(@coverageerr,startph,[],[],options,Theta,Phi,-bwim);

% if resnorm2<resnorm
%     fitp=fitp2;
% end
if nargin>2
imfit=startp(4)*(coverage_sphere(fitp(1),fitp(2),fitp(3),Theta,Phi)+1);
% im=coverage(-pi/4,pi/4,pi/16,Theta,Phi);
imcombine=zeros(size(imstart,2),size(imstart,1),3);
imcombine(:,:,1)=imin'/max(imin(:))*3;
imcombine(:,:,3)=imfit';
imcombine(:,:,2)=imcombine(:,:,1);
% imcombine(:,:,2)=imstart';
imagesc(ax,pr,tr,imcombine);
title(ax,['Image: Theta=' num2str(180-fitp(3)/pi*180,'%3.0f') '°'])
% figure(81);imagesc(tr,pr,imfit);
end

end

function err=coverageerr(fitp,Theta,Phi,imin)
% fitp
% fitp(3)=pi/4;
imtest=coverage_sphere(fitp(1),fitp(2),fitp(3),Theta,Phi);
% err=(imtest(:)+1)*fitp(4)+fitp(5)-imin(:);
err=(imtest(:)+1)*fitp(4)-imin(:);
if 0
imcombine=zeros(size(imin,1),size(imin,2),3);
imcombine(:,:,1)=(imin+1)/2;
imcombine(:,:,2)=(imtest+1)/2;
imcombine(:,:,3)=imstart;
figure(77);imagesc(imcombine);title(fitp)
drawnow
end
end