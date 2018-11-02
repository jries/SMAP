function [fitout,fout,fstart]=convoluteshapefits(xv,profile,sigma,structure)
profile=double(profile);
xv=double(xv);
sigma=double(sigma);

[m,indm]=max(profile);

s1=xv(find(profile>m/2,1,'first'));
s2=xv(find(profile>m/2,1,'last'));
x0s=(s1+s2)/2;
oldopts=optimset('lsqcurvefit');
newopts=optimset(oldopts,'MaxFunEvals',150,'MaxIter',150,'TolFun',1e-3,'Display','off');%,'Algorithm','levenberg-marquardt');


switch structure
    
    case 1 %stepfunction
        start=[m,s1,s2-s1,0,sigma]
        fitout=lsqcurvefit(@fitstep,double(start),double(xv),double((profile)),[],[],newopts);
        fout=fitstep(fitout,xv);
        fstart=fitstep(start,xv);
        
   case 2 %disk

%         a0=fitdisk([1,x0s,(s2-s1)/2,0],x0s,sigma);
        start=[m,x0s,(s2-s1)/2,mean(profile(1:3)),sigma];
        fitout=lsqcurvefit(@fitdisk,double(start),double(xv),double((profile)),[],[],newopts);
        
        fout=fitdisk(fitout,xv);  
        fstart=fitdisk(start,xv);
%           fout=fitdisk(start,xv,sigma); 
   case 3 %ring
        start=[m,x0s,(s2-s1)/2,mean(profile(1:3)),sigma];
        fitout=lsqcurvefit(@fitring,double(start),double(xv),double((profile)),[],[],newopts);
        fout=fitring(fitout,xv);
         fstart=fitring(start,xv);
   case 4 %two structure: distance. Gaussian Approximation
       start=[m,m,s1,s2-s1,sigma,0];
        fitout=lsqcurvefit(@fit2gauss,double(start),double(xv),double((profile)),[],[],newopts);
        fout=fit2gauss(fitout,xv);
         fstart=fit2gauss(start,xv);

end

function out=fit2gauss(p,xv)
%p=[A,x0,R], 
A1=p(1);
A2=p(2);
x0=p(3);
d=p(4);
s=p(5);
off=p(6);
out=A1*exp(-(xv-x0).^2/2/s^2)+A2*exp(-(xv-x0-d).^2/2/s^2)+off;


function out=fitstep(p,xv)
%p=[A,x0,L], 
A=p(1);
x0=p(2);
L=p(3);
off=p(4);
sigma=p(5);
out=(quadv(@(x)convstep(x+x0,xv,sigma),0,L));
out=A*out/max(out)+off;

function out=fitring(p,xv,sigma)
%p=[A,x0,R], 

A=p(1);
x0=p(2);
R=max(0.00001,p(3));
off=p(4);
sigma=p(5);
out=quadv(@(x)convring(x,xv,sigma,R,x0),-R+x0,R+x0);
out=A*out/max(out)+off;

function out=fitdisk(p,xv,sigma)
%p=[A,x0,R], 
A=p(1);
x0=p(2);
R=max(0.00001,p(3));
off=p(4);
sigma=p(5);
out=quadv(@(x)convdisk(x,xv,sigma,R,x0),-R+x0,R+x0,1e-4);
out=out/max(out)*A+off;


function Y = convstep(xs,x,pp,x0)
Y=profile(x-xs,pp);
function Y = convring(xs,x,pp,ps,x0)
Y=profile(x-xs,pp).*structurering(xs-x0,ps);

function Y = convdisk(xs,x,pp,ps,x0)
Y=profile(x-xs,pp).*structuredisk(xs-x0,ps);

function out=profile(x,pp)
out=exp(-x.^2/2/pp(1)^2);

function out=structuredisk(x,ps)
out=sqrt(ps(1)^2-x.^2);
if imag(out)
    out=-1000+0*x;
end

function out=structurering(x,ps)
out=ps(1)/sqrt(ps(1)^2-x.^2);
if imag(out)
    out=-1000+0*x;
end