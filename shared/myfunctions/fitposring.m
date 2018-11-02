function [x0,y0,R0]=fitposring(x,y,R,startpar)
if nargin<3||isempty(R) %fit also R
    fh=@freering;
    if nargin<4    
        xs=mean(x);ys=mean(y);
        rs=sqrt(std(x).^2+std(y).^2);
        startpar=[xs,ys,rs];
    end
    fitp=implicitfit(fh,startpar,x,y,0);
    x0=fitp(1);y0=fitp(2);R0=fitp(3);
else 
    fh=@fixring;
    if nargin<4    
        xs=mean(x);ys=mean(y);
        startpar=[xs,ys];
    end
    fitp=implicitfit(fh,startpar,x,y,R);
    x0=fitp(1);y0=fitp(2);
    R0=R;
end


function err=fixring(par,x,y,R,d1,d2,d3)
err=sqrt(((x-par(1)).^2)+((y-par(2)).^2))-R;

function err=freering(par,x,y,z,d1,d2,d3)
R=par(3);
err=sqrt(((x-par(1)).^2)+((y-par(2)).^2))-R;