function out=imgringfit(img)
% img is quare

roisize=size(img,1);

fit_sigma=true; %also fit the Gaussian blurr
sigma=2; %magnitude of gaussian blurr, or startparameter therefore
x0=0;y0=0; %startparameters
r0=5; %startradius
startdr=1; %start ring thickness
exponent=1; %dont re normalize image.
% maxradius=100;
n=-floor(roisize/2):-floor(roisize/2)+roisize-1;
 [X,Y]=meshgrid(n);

 
 
 img=(img).^(1/exponent);
  startp=[mean(img(:)), x0,y0,r0+startdr,startdr];
  
 
  lb=[0 -inf -inf 0 0];
  ub=[inf inf inf roisize*2 roisize*2];
 
  if fit_sigma
   startp(end+1)=sigma;
   lb(end+1)=sigma/1.5;
   ub(end+1)=sigma*3;
  end
   
 fitpim=lsqcurvefit(@gaussring,startp,X,img,lb,ub,[],Y,sigma,exponent);
 imfit=gaussring(fitpim,X,Y,sigma,exponent);
 imstart=gaussring(startp,X,Y,sigma,exponent);
 
 
 figure(88);
 subplot(2,1,1)

 imagesc(horzcat(img,horzcat(imstart,imfit)));
 subplot(2,1,2)
  hold off
 imagesc(n,n,img)
hold on
   circle(fitpim(2),fitpim(3),fitpim(4),'EdgeColor','m')
    circle(fitpim(2),fitpim(3),fitpim(4)-fitpim(5),'EdgeColor','m')
   axis equal
   out.x0=fitpim(2);out.y0=fitpim(3);
   out.r1=fitpim(4);out.dr1=fitpim(5);
   
   if length(fitpim)>5
       out.sigma1=fitpim(6);
   else
       out.sigma1=sigma;
   end
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
R=sqrt((X-x).^2+(Y-y).^2);
im=(a*(erf((r-R)/sigma2)-erf(((r-dr)-R)/sigma2)));
 im=(im).^(1/exponent);
end

function circle(x,y,r,varargin)
if r>0
posf=[x-r,y-r,2*r,2*r];
rectangle('Position',posf,'Curvature',[1 1],varargin{:})
end
end