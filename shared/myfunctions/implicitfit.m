function [fitp,resnorm]=implicitfit(fh,startp,x,y,z,lb,ub,coarse)
if nargin<6
lb=[];
ub=[];
end
if nargin<8
coarse=false;
end
wx=0;wy=0;wz=0;
 opt=optimset('lsqnonlin');
 opt.Display='off';
 if coarse
 opt.TolFun=3e-3;
 opt.MaxFunEvals=50;
 end
[fitp,residual,~,exitflag]=lsqnonlin(@callfit,double(startp),lb,ub,opt,fh,double(x),double(y),double(z),double(wx),double(wy),double(wz));

resnorm=sqrt(sum(residual.^2));


function err=callfit(par,fh,x,y,z,wx,wy,wz)
%     err=sqrt(abs(fh(par,x,y,z,wx,wy,wz)));
    err=((fh(par,x,y,z,wx,wy,wz)));
err=err./sqrt(abs(err));