fn='/Volumes/LaCie/otherData/deepSMLM/live/beadstacks561nm2DnoAirbubble/nomirror_beadstacks561nm2DnoAirbubble_3dcal.mat';
l=load(fn);

imgstack=l.SXY.PSF{1};
imgc=double(imgstack(2:end-1,2:end-1,:));

mp=(size(imgc,1)-1)/2;

[X,Y]=meshgrid(-mp:mp);
[Xn,Yn]=meshgrid(-mp-1:mp+1);
g1=fittype( @(a, b, c, d, x, y) a*exp(-((x-b).^2+(y-c).^2)/d^2), 'independent', {'x', 'y'},...
     'dependent', 'z' );
g2=fittype( @(a, b, c, d,e,f, x, y) a*exp(-((x-b).^2+(y-c).^2)/d^2)+e*exp(-((x-b).^2+(y-c).^2)/f^2), 'independent', {'x', 'y'},...
     'dependent', 'z' );
g3=fittype( @(a, b, c, d,e,f, g,h,x, y) a*exp(-((x-b).^2+(y-c).^2)/d^2)+e*exp(-((x-b).^2+(y-c).^2)/f^2)+g*exp(-((x-b).^2+(y-c).^2)/h^2), 'independent', {'x', 'y'},...
     'dependent', 'z' );
g4=fittype( @(a, b, c, d,e,f, g,h,m,n,x, y) a*exp(-((x-b).^2+(y-c).^2)/d^2)+e*exp(-((x-b).^2+(y-c).^2)/f^2)+g*exp(-((x-b).^2+(y-c).^2)/h^2)+m*exp(-((x-b).^2+(y-c).^2)/n^2), 'independent', {'x', 'y'},...
     'dependent', 'z' );
 
 
 mps=ceil(size(imgc,3)/2);
 imh=imgc(:,:,mps);
 maxim=max(imh(:));
 sp2_0=[maxim, 0,0,1,0,3];
 sp3_0=[maxim, 0,0,1,0,3,0,5];
     
  sp4_0=[maxim, 0,0,1,0,3,0,5,0,15];
  
  sp1=[maxim, 0,0,3];
  
 l2=[-inf -2 -2 0.5 -inf 2];
 u2=[inf 2 2 3 inf 30];
 
  l3=[-inf -2 -2 1.5 -inf 4 -inf 0.7];
 u3=[inf 2 2 4 inf 30 inf 1.5];
 
 
 sp2=sp2_0;
 sp3=sp3_0;
 
%  sp4=[maxim/2, 0,0,10,-maxim/5,6,maxim/2,3,-maxim/5,1.5];
%  
%  for k=161:-1:1
%     imh=imgc(:,:,k);
% 
% %     fitp1=fit([X(:),Y(:)],imh(:),g1,'StartPoint',sp1);
% %     s1(k)=fitp1.d;
% %     
% %      sp2=[fitp1.a, fitp1.b,fitp1.c,fitp1.d,0,fitp1.d*3];
% %      sp3=[fitp1.a, fitp1.b,fitp1.c,fitp1.d,0,fitp1.d*3,0,1.5];
% %      sp4=[fitp1.a, fitp1.b,fitp1.c,fitp1.d,0,fitp1.d*2,0,1.5,0,fitp1.d/2];
% %  
% %         fitp2=fit([X(:),Y(:)],imh(:),g2,'StartPoint',sp2);
% %     sp2=coeffvalues(fitp2);
% %     fitp3=fit([X(:),Y(:)],imh(:),g3,'StartPoint',sp3);
% %     sp3=coeffvalues(fitp3);
%         fitp4=fit([X(:),Y(:)],imh(:),g4,'StartPoint',sp4);
%     sp4=coeffvalues(fitp4);
% %     subplot(2,1,1);imagesc()
% %     implot2(:,:,k)=[imh, fitp2(X,Y);imh-fitp2(X,Y),imh*0];
% %  implot3(:,:,k)=[imh, fitp3(X,Y);imh-fitp3(X,Y),imh*0]; 
%  implot4(:,:,k)=[imh, fitp4(X,Y);imh-fitp4(X,Y),imh*0];
%     
%  end
%  
%  
%  imx(vertcat(implot2,implot3,implot4))

 %
 sp4=[maxim/2, 0,0,10,-maxim/5,6,maxim/2,3,-maxim/5,1.5];
 mps=161;
 clear psf4
for k=mps:-1:1
    imh=imgc(:,:,k);

%     fitp2=fit([X(:),Y(:)],imh(:),g2,'StartPoint',sp2,'Lower',l2,'Upper',u2);
%     sp2=coeffvalues(fitp2);
%     fitp3=fit([X(:),Y(:)],imh(:),g3,'StartPoint',sp3,'Lower',l3,'Upper',u3);
%     sp3=coeffvalues(fitp3);
% %     subplot(2,1,1);imagesc()
%     implot2(:,:,k)=[imh, fitp2(X,Y);imh-fitp2(X,Y),imh*0];
%  implot3(:,:,k)=[imh, fitp3(X,Y);imh-fitp3(X,Y),imh*0];
%  
         fitp4=fit([X(:),Y(:)],imh(:),g4,'StartPoint',sp4);
    sp4=coeffvalues(fitp4);
     implot4(:,:,k)=[imh, fitp4(X,Y);imh-fitp4(X,Y),imh*0];
     fitpc=fitp4;
     fitpc.b=0;fitpc.c=0;
     psf4(:,:,k)=fitpc(Xn,Yn);
end
 sp2=sp2_0;
 sp3=sp3_0;
  sp4=[maxim/2, 0,0,10,-maxim/5,6,maxim/2,3,-maxim/5,1.5];
for k=mps:size(imgc,3)
    imh=imgc(:,:,k);

%     fitp2=fit([X(:),Y(:)],imh(:),g2,'StartPoint',sp2);
%     sp2=coeffvalues(fitp2);
%     fitp3=fit([X(:),Y(:)],imh(:),g3,'StartPoint',sp3);
%     sp3=coeffvalues(fitp3);
% %     subplot(2,1,1);imagesc()
%     implot2(:,:,k)=[imh, fitp2(X,Y);imh-fitp2(X,Y),imh*0];
%  implot3(:,:,k)=[imh, fitp3(X,Y);imh-fitp3(X,Y),imh*0];
          fitp4=fit([X(:),Y(:)],imh(:),g4,'StartPoint',sp4);
    sp4=coeffvalues(fitp4);
     implot4(:,:,k)=[imh, fitp4(X,Y);imh-fitp4(X,Y),imh*0];
     fitpc=fitp4;
     fitpc.b=0;fitpc.c=0;
     psf4(:,:,k)=fitpc(Xn,Yn);
end

imx(implot4)
imx(psf4)


psfold=splinePSF;
psfold.loadmodel(fn);

zp=(-psfold.modelpar.z0+1:psfold.modelpar.z0-1)*psfold.modelpar.dz;
coordh(:,3)=-zp;
imagex=psfold.PSF(coordh,31);
figure(33);hold off
plot(squeeze(imagex(16,16,:)));
hold on
plot(squeeze(psf4(18,18,:)))

figure(34)
hold off
plot(squeeze(sum(sum(psf4n,1),2)))
hold on
plot(squeeze(sum(sum(imagex,1),2)))

psf4n=psf4/sum(sum(psf4(:,:,round(end/2))));

coeff = Spline3D_interp(psf4n(:,:,2:end-1));

lout=l;
lout.SXY.PSF{1}=psf4;
lout.SXY.cspline.coeff{1}=coeff;
outf=strrep(fn,'_3dcal','_4Gauss_3dcal');
save(outf,'-struct','lout')

%%
psf=splinePSF;
psf.loadmodel(outf);
z=(-700:1:700)';
profilep=psf.PSF([0*z,0*z,z],1);
figure(80);hold off;plot(z,squeeze(profilep))
hold on
profileo=psfold.PSF([0*z,0*z,z],1);
plot(z,squeeze(profileo))
