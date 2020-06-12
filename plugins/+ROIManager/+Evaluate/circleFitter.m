classdef circleFitter<interfaces.SEEvaluationProcessor
%     Fits circular model to localizations in the ROI to determine the
%     radius and center coordinates of circular structures.
    properties
        boundary
    end
    methods
        function obj=circleFitter(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function out=run(obj, p)
            layers=find(obj.getPar('sr_layerson'));
            for layer=layers
                out=runintern(obj,p,layer);
            end
        end
     
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end

end


function out=runintern(obj,p,layer)
p=obj.getLayerParameters(layer);
if p.sqrtfit
    p.exponent=2;
else
    p.exponent=1;
end
roisize=p.se_siteroi/2;

locs=obj.getLocs({'xnm','ynm','locprecnm'},'layer',layer,'size',roisize);       
xm=locs.xnm-obj.site.pos(1); %coordinates relative to site position
ym=locs.ynm-obj.site.pos(2); 
ho=obj.setoutput(['circle' num2str(layer)]);
 % circfit
circfit=ringfit(p,xm,ym,ho);
%Ng and Nu
locsg=obj.getLocs({'xnmr','ynmr'},'layer',layer,'size',roisize,'grouping','grouped');
locsu=obj.getLocs({'xnmr','ynmr'},'layer',layer,'size',roisize,'grouping','ungrouped');
circfactor=1.5;
circfit.Ncirc1=sum((locs.xnmr-circfit.x0).^2+(locs.ynmr-circfit.y0).^2<circfit.r1.^2*circfactor^2);
circfit.Ncircg1=sum((locsg.xnmr-circfit.x0).^2+(locsg.ynmr-circfit.y0).^2<circfit.r1.^2*circfactor^2);
circfit.Ncircu1=sum((locsu.xnmr-circfit.x0).^2+(locsu.ynmr-circfit.y0).^2<circfit.r1.^2*circfactor^2);


xf=xm-circfit.x0;yf=ym-circfit.y0;


%imfit
hf=obj.setoutput(['image' num2str(layer)]);
imfit=imgfit(p,xm,ym,locs,circfit,hf,layer);
imfit.Ncirc1=sum((locs.xnmr-imfit.x0).^2+(locs.ynmr-imfit.y0).^2<imfit.r1.^2*circfactor^2);
imfit.Ncircg1=sum((locsg.xnmr-imfit.x0).^2+(locsg.ynmr-imfit.y0).^2<imfit.r1.^2*circfactor^2);
imfit.Ncircu1=sum((locsu.xnmr-imfit.x0).^2+(locsu.ynmr-imfit.y0).^2<imfit.r1.^2*circfactor^2);
xf=xm-imfit.x0;yf=ym-imfit.y0;


imfit.asymmetry=getasymmetry(p,xf,yf);
out.imfit=imfit;out.circfit=circfit;

%make larger image
locs=obj.getLocs({'xnm','ynm','xnmr','ynmr','locprecnm'},'layer',layer,'size',[p.se_sitefov p.se_sitefov]);

xm=locs.xnm-obj.site.pos(1);
ym=locs.ynm-obj.site.pos(2);

xc=xm-imfit.x0;
yc=ym-imfit.y0;
s=locs.locprecnm/2;
out.imfit.image=makeimage(p,xc,yc,s);

end

function im=makeimage(p,xm,ym,s)
roisize=p.se_siteroi/2;
pixels=p.se_sitepixelsize;
 range=[-roisize roisize];
 posf.x=xm;posf.y=ym;posf.s=s;
im=double(gaussrender(posf,range, range, pixels, pixels));
end

function out=getasymmetry(p,x,y)
I(1,1)=sum(x.^2);
I(2,2)=sum(y.^2);
I(1,2)=-sum(x.*y);
I(2,1)=I(1,2);
I=I/length(x);

e=(eig(I));

out=sort(sqrt(e));
end

function out=imgfit(p,xm,ym,locs,circfit,hf2,layer)
roisize=p.se_siteroi/2;
 posf.x=xm;posf.y=ym;
 posf.s=locs.locprecnm*p.gaussfac*p.gaussfac_imfit(min(layer,length(p.gaussfac_imfit)));
 ming=max(p.mingaussnm,p.mingausspix*p.se_sitepixelsize) 	;
 posf.s(posf.s<ming)=ming;
 range=[-roisize roisize];
 pixels=p.se_sitepixelsize;
 img=double(gaussrender(posf,range, range, pixels, pixels));

 ming2=max(ming,median(posf.s));
% max( max(gaussrender(struct('x',0,'y',0,'s',ming),[-10 10], [-10 10], pixels, pixels)))
 
% THRESHOLD for the fit, any intensity values above co are saturated for
% the fit
% co=1/(2*pi*ming2^2/pixels^2)*4; % uses minimum gaussian setting
% co=myquantilefast(img,.95)*3; % uses 95th quantile of intensities
% img(img>co)=co;
%
 
%  qim=myquantile(img(:),.95);
%  img(img>qim)=qim;
startdr=30;
maxradius=100;
 [X,Y]=meshgrid(-roisize+pixels/2:pixels:roisize);
  startpc=circfit;
 if startpc.r1>maxradius
     startpc.x0=0;startpc.y0=0;startpc.r1=roisize/4;
 end
 startdr=startpc.r1*.75;
 

 sigma=double(median(locs.locprecnm))*sqrt(1+(p.gaussfac_imfit(1)*p.layer1_.gaussfac)^2);%*sqrt(p.exponent);
 imgf=(img).^(1/p.exponent);
  startp=[mean(img(:)), startpc.x0,startpc.y0,startpc.r1+startdr/2,startdr];
  
 
  lb=[0 -inf -inf 0 0];
  ub=[inf inf inf roisize*2 roisize*2];
 
  if p.fit_sigma
   startp(end+1)=sigma;
   lb(end+1)=sigma/1.5;
   ub(end+1)=sigma*3;
  end
   
 fitpim=lsqcurvefit(@gaussring,startp,X,imgf,lb,ub,[],Y,sigma,p.exponent);
 imfit=gaussring(fitpim,X,Y,sigma,p.exponent);
 imstart=gaussring(startp,X,Y,sigma,p.exponent);
%   imagesc(vertcat(img.^2,imfit.^2,img.^2-imfit.^2),'Parent',hf)
%     imagesc(vertcat(horzcat(img,imfit),horzcat(imstart,img-imfit)),'Parent',hf)
range2=[-roisize 3*roisize];

% I0=2*pi*sigma^2/pixels^2;

imres=imgf.^p.exponent-imfit.^p.exponent;
imres=imres/max(imres(:));
imgf=imgf/max(imgf(:));
imfit=imfit/max(imfit(:));
imstart=imstart/max(imstart(:));
%    imagesc(range2,range2,vertcat(horzcat(imgf.^p.exponent,imfit.^p.exponent),horzcat(imstart.^p.exponent,imres)),'Parent',hf)
   imagesc(range,range,imgf.^p.exponent,'Parent',hf2)

   hf2.NextPlot='add';
   circle(circfit.x0,circfit.y0,circfit.r1,'Parent',hf2);
   circle(fitpim(2),fitpim(3),fitpim(4),'EdgeColor','m','Parent',hf2)
%    if fitpim(5)<fitpim(4)
    circle(fitpim(2),fitpim(3),fitpim(4)-fitpim(5),'EdgeColor','m','Parent',hf2)
%    end
%    circle(fitpim(2),fitpim(3),fitpim(4),'EdgeColor','r','Parent',hf2)
   hf2.NextPlot='replace';
   xlabel(hf2,'x (nm)');
   ylabel(hf2,'y (nm)');
   
   fm='%2.1f';
   if p.fit_sigma
       tsig=[' nm, \sigma = ' num2str(fitpim(6),fm) ' nm'];
       out.sigma1=fitpim(6);
   else
       tsig=[' nm, fix \sigma = ' num2str(sigma,fm) ' nm'];
       out.sigma1=sigma;
   end
   
   title(['R_o = ' num2str(fitpim(4),fm) ' nm, d_r = ' num2str(fitpim(5),fm) tsig]);
   
   out.x0=fitpim(2);out.y0=fitpim(3);
   out.r1=fitpim(4);out.dr1=fitpim(5);
   out.dr_ro1=out.dr1/(out.r1);
   
%    out.image=img.^2;
  %segmentation:
  h=fspecial('gaussian',7,3);
  imh=filter2(h,img);
  imbw=imh>max(h(:))*1.5;
  imbws=bwareafilt(imbw,1);
  try
  out.areastat=regionprops(imbws,'MajorAxisLength','MinorAxisLength','ConvexArea','Area','Eccentricity');
  out.areastat.fractionCircle=out.areastat.Area/(pi*(out.areastat.MajorAxisLength/2)^2);
  out.areastat.fractionCircleconv=out.areastat.ConvexArea/(pi*(out.areastat.MajorAxisLength/2)^2);
  catch err
      out.areastat(1).MajorAxisLength=0;
      out.areastat(1).MinorAxisLength=0;
      out.areastat(1).ConvexArea=0;
      out.areastat(1).Area=0;
      out.areastat(1).Eccentricity=0;
      out.areastat(1).fractionCircle=0;
      out.areastat(1).fractionCircleconv=0;
  end

end
function out=imgfit2(p,xm1,ym1,xm2,ym2,locs1,locs2,circfit,hf,hf2)
roisize=p.se_siteroi/2;
 posf1.x=xm1;posf1.y=ym1;posf1.s=locs1.locprecnm*p.layer1_.gaussfac*p.gaussfac_imfit(1);
 posf2.x=xm2;posf2.y=ym2;posf2.s=locs2.locprecnm*p.layer2_.gaussfac*p.gaussfac_imfit(2);
 range=[-roisize roisize];
 pixels=p.se_sitepixelsize;
 disp('render')
 img1=double(gaussrender(posf1,range, range, pixels, pixels));
 img2=double(gaussrender(posf2,range, range, pixels, pixels));
%  qim=myquantile(img(:),.95);
%  img(img>qim)=qim;
startdr=30;
maxradius=75;%for startp parameter only
 [X,Y]=meshgrid(-roisize+pixels/2:pixels:roisize);
  startpc=circfit;
 if startpc.r1>maxradius||startpc.r2>maxradius
     startpc.x0=0;startpc.y0=0;startpc.r1=roisize/4;startpc.r2=roisize/4;
 end
 

 sigma1=double(median(locs1.locprecnm))*sqrt(1+p.gaussfac_imfit(1)^2*p.layer1_.gaussfac^2)*sqrt(p.exponent);
 sigma2=double(median(locs2.locprecnm))*sqrt(1+p.gaussfac_imfit(2)^2*p.layer2_.gaussfac^2)*sqrt(p.exponent);
 
 img1=(img1).^(1/p.exponent);
 img2=(img2).^(1/p.exponent); 
%   disp('imgfit2 fit') 
opt=optimset('lsqnonlin');
opt.MaxIter=100;
opt.TolX=1e-4;
% opt.Display='off';
opt.TolFun=1e-4;

astart=[max(img1(:))/2,max(img2(:))/2];
astart=[mean(img1(:)),mean(img2(:))];
 startp=[astart(1),astart(2) startpc.x0,startpc.y0,startpc.r1+startdr/2,startpc.r2+startdr/2,startdr,startdr];
 fitpim=lsqnonlin(@gaussringerr,startp,[0*startp(1)/5 0*startp(2)/5 -inf -inf 0 0 sigma1/4 sigma2/4],...
     [inf inf inf inf roisize*2 roisize*2 roisize*2 roisize*2],opt,img1,img2,X,Y,sigma1,sigma2,p.exponent);
%  disp('imgfit2 fit done') 
 p1=fitpim([1 3 4 5 7]);
p2=fitpim([2 3 4 6 8]);
ps1=startp([1 3 4 5 7]);
ps2=startp([2 3 4 6 8]);
imfit1=gaussring(p1,X,Y,sigma1,p.exponent);
imfit2=gaussring(p2,X,Y,sigma2,p.exponent);
imfit=vertcat(imfit1/max(imfit1(:)),imfit2/max(imfit2(:))).^p.exponent;

imstart1=gaussring(ps1,X,Y,sigma1,p.exponent);
imstart2=gaussring(ps2,X,Y,sigma2,p.exponent);
imstart=vertcat(imstart1/max(imstart1(:)),imstart2/max(imstart2(:))).^p.exponent;

images=vertcat(img1/max(img1(:)),img2/max(img2(:))).^p.exponent;
imagesside=horzcat(img1,img2).^p.exponent;
%   imagesc(vertcat(img.^2,imfit.^2,img.^2-imfit.^2),'Parent',hf)
%     imagesc(vertcat(horzcat(img,imfit),horzcat(imstart,img-imfit)),'Parent',hf)
range3=[-roisize 5*roisize];
range2=[-roisize 3*roisize];
   imagesc(range3,range2,horzcat(images,imfit,imstart),'Parent',hf)
   imagesc(range2,range,imagesside,'Parent',hf2)

   hf2.NextPlot='add';
   circle(circfit.x0,circfit.y0,circfit.r1,'Parent',hf2);
   circle(circfit.x0+2*roisize,circfit.y0,circfit.r2,'Parent',hf2);
   circle(p1(2),p1(3),p1(4)-p1(5),'EdgeColor','m','Parent',hf2)
   circle(p2(2)+2*roisize,p2(3),p2(4)-p2(5),'EdgeColor','m','Parent',hf2)
%    if fitpim(5)<fitpim(4)
%     circle(fitpim(2),fitpim(3),fitpim(4)-fitpim(5),'EdgeColor','m','Parent',hf2)
%    end
   circle(p1(2),p1(3),p1(4),'EdgeColor','r','Parent',hf2)
   circle(p2(2)+2*roisize,p2(3),p2(4),'EdgeColor','r','Parent',hf2)
   hf2.NextPlot='replace';
   
   out.x0=p1(2);out.y0=p1(3);
   out.r1=p1(4);out.dr1=p1(5);
   out.r2=p2(4);out.dr2=p2(5);
   out.dr_ro1=out.dr1/(out.r1);
   out.dr_ro2=out.dr2/(out.r2);
   out.sigma1=sigma1;
   out.sigma2=sigma2;
   
   %segment:
   
   
%    out.image=img.^2;
%   out
%  fitpim(1:2)
%  startp(1:2)
end

function out=rthetaprofiles(p,xf,yf,hr,ht)
%  figure(88);plot(xf,yf,'+')
roisize=p.se_siteroi/2;
[theta,rho]=cart2pol(xf,yf);

 dr=5;
%  rn=dr/2:dr:roisize;
rn=0:dr:roisize;
 histr=histcounts(rho,rn);
 
 norm=rn(2:end).^2-rn(1:end-1).^2;
%  rdensity=histr./(rn+dr/2);
rdensity=histr./norm;
 plot(rn(1:end-1)+dr/2,rdensity,'Parent',hr)
 out.rdensity=rdensity;
 out.rn=rn;

 dtheta=pi/16;
 thetan=-pi+dtheta/2:dtheta:pi;
 histtheta=hist(theta,thetan);
  histtheta2=hist(theta+pi,thetan+pi);
  
 tac=myxcorr(histtheta,histtheta);
 tac=tac+myxcorr(histtheta2,histtheta2);
 plot(thetan+pi,tac,'Parent',ht)
 
 %look at variance along theta
 dtheta=2*pi/8; %coarse
 thetanh=-pi+dtheta/2:dtheta:pi;
 histtheta=hist(theta,thetanh);
 out.theta_varnorm=var(histtheta)/mean(histtheta);
 out.xcorramp=tac(1)/(mean(histtheta)+mean(histtheta2))^2;
%  figure(88);plot(thetanh,histtheta); title([out.theta_varnorm out.xcorramp]);
 
 out.thetaAC=tac;
 out.thetan=thetan+pi;
 out.histtheta=histtheta;
end

function out=ringfit2(p,x1,y1,x2,y2,ho)
roisize=p.se_siteroi/2;
start_radius=roisize/2;
startp=[0 0 start_radius start_radius];
lb=[-inf -inf 0 0];
ub=[inf inf roisize*2 roisize*2];
[fitp,residual,~,exitflag]=lsqnonlin(@fit2circle,double(startp),lb,ub,[],double(x1),double(y1),double(x2),double(y2));
plot(x1,y1,'x','Parent',ho)
ho.NextPlot='add';
plot(x2,y2,'x','Parent',ho)
circle(fitp(1),fitp(2),fitp(3),'Parent',ho);
circle(fitp(1),fitp(2),fitp(4),'Parent',ho);
% pos = [fitp(2)-fitp(1) fitp(3)-fitp(1) 2*fitp(1) 2*fitp(1)];
% rectangle('Position',pos,'Curvature',[1 1],'Parent',ho)
ho.NextPlot='replace';
axis(ho,'ij','equal');
out.x0=fitp(1);
out.y0=fitp(2);
out.r1=fitp(3);
out.r2=fitp(4);
end
function error=fit2circle(par,x1,y1,x2,y2)
error=vertcat(sqrt(((x1-par(1)).^2)+((y1-par(2)).^2))-par(3),sqrt(((x2-par(1)).^2)+((y2-par(2)).^2))-par(4));
end
function out=ringfit(p,xm,ym,ho)
roisize=p.se_siteroi/2;
fh=@circle_implicit;
start_radius=roisize/2;
startp=[start_radius,mean(xm),mean(ym)];

lb=[0 -2 -2 ]*1000;
ub=[2 2 2 ]*1000;
% whos xm
[fitp,r]=implicitfit(fh,startp,xm,ym,0*xm,lb,ub);
% [xy, rad]=fitcircle([xm,ym]');
% fitp2=[xy(1) xy(2) rad];

plot(xm,ym,'x','Parent',ho)
ho.NextPlot='add';
% circle(fitp2(1),fitp2(2),fitp2(3),'EdgeColor',[1 0 0],'Parent',ho);
circle(fitp(2),fitp(3),fitp(1),'Parent',ho);
% pos = [fitp(2)-fitp(1) fitp(3)-fitp(1) 2*fitp(1) 2*fitp(1)];
% rectangle('Position',pos,'Curvature',[1 1],'Parent',ho)
ho.NextPlot='replace';
axis(ho,'ij','equal');
xlabel(ho,'x (nm)');
ylabel(ho,'y (nm)');

 circfactor=1.5;
out.r1=fitp(1);
out.x0=fitp(2);
out.y0=fitp(3);
out.Ncirc=sum((xm-fitp(2)).^2+(ym-fitp(3)).^2<fitp(1).^2*circfactor^2);
title(['N=' num2str(out.Ncirc) ', r=' num2str(out.r1)],'Parent',ho)
end

function err=gaussringerr(par,img1,img2,X,Y,s1,s2,exponent)
%a1 a2 x y r1 r2 dr1 dr2
p1=par([1 3 4 5 7]);
p2=par([2 3 4 6 8]);

g1=gaussring(p1,X,Y,s1,exponent);
g2=gaussring(p2,X,Y,s2,exponent);
% err1=(g1-img1).*(sqrt(img1)+mean(img1(:)));
% err2=(g2-img2).*(sqrt(img2)+mean(img2(:)));
err1=(g1-img1)*sum(img1(:));
err2=(g2-img2)*sum(img2(:));
err=vertcat(err1(:),err2(:));
end

function im=gaussring(par,X,Y,sigma,exponent)

a=par(1);
x=par(2);
y=par(3);
r=par(4);
dr=par(5);
if length(par)>5
    sigma=par(6);
end
sigma2=sigma*sqrt(2);
% dr=10;

R=sqrt((X-x).^2+(Y-y).^2);
im=(a*(erf((r-R)/sigma2)-erf(((r-dr)-R)/sigma2)));
% im=(a*(erf((R+dr-r)/sigma)-erf(((R-dr)-r)/sigma)));
% min(im(:))
% max(im(:))
 im=(im).^(1/exponent);
end

function pard=guidef(obj)
% p(1).value=1; p(1).on={};p(1).off={'sqrtfit','gaussfac_imfitt','gaussfac_imfit','fit_sigma'};
% p(2).value=2; p(2).off={};p(2).on=p(1).off;
% pard.fitmode.object=struct('Style','popupmenu','String',{{'Circle','Ring'}},'Value',2,'Callback',{{@obj.switchvisible,p}});
% pard.fitmode.position=[1,1];
% pard.fitmode.Width=2;

pard.sqrtfit.object=struct('Style','checkbox','String','fit sqrt(img)','Value',1);
pard.sqrtfit.position=[2,1];
pard.sqrtfit.Width=2;

pard.gaussfac_imfitt.object=struct('Style','text','String','gaussfac for imfit');
pard.gaussfac_imfitt.position=[3,1];
pard.gaussfac_imfitt.Width=2;

pard.gaussfac_imfit.object=struct('Style','edit','String','1');
pard.gaussfac_imfit.position=[3,3];
pard.gaussfac_imfit.Width=1;

pard.fit_sigma.object=struct('Style','checkbox','String','Fit sigma of ring','Value',0);
pard.fit_sigma.position=[4,1];
pard.fit_sigma.Width=2;

pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi'};
pard.plugininfo.type='ROI_Evaluate';
pard.plugininfo.description='Fits circular model to localizations in the ROI to determine the radius and center coordinates of circular structures.';
end
