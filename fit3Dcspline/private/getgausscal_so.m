function [gauss]=getgausscal_so(curves,p)
        p.ax=p.ax_sxsy;
        gauss.Sx2_Sy2=cal_Sx2_Sy2(curves,p);
   
        p.ax=p.ax_z;
        gauss.fitzpar=cal_fitzpar(curves,p);
    drawnow
end



function sxp=cal_Sx2_Sy2(b,p)


z=vertcat(b(:).z);
Sx=vertcat(b(:).sx);
Sy=vertcat(b(:).sy);
hold off
sxp=fitsx2sy2_so(Sx,Sy,z,p.gaussrange,p.ax);
xlim(p.ax,[p.gaussrange])
end

function fitzpar=cal_fitzpar(b,p)
 z=vertcat(b(:).z);
Sx=vertcat(b(:).sx);
Sy=vertcat(b(:).sy);
zrange=[p.gaussrange(1) p.gaussrange(2)];
fitzpar=getzfitpar(Sx,Sy,z,zrange,0,true,p.ax);
end

function out=fitsx2sy2_so(sx,sy,z,zrange,ax)
indf=z>zrange(1)&z<zrange(2);

ds=sx.^2-sy.^2;
q=myquantile(ds(indf),[0.05 0.95])+[-2 2];

inds=ds>q(1)&ds<q(2);
inds=indf&inds;

if sum(indf)>5
warning('off','curvefit:fit:iterationLimitReached');
fitpsx=fit(sx(inds).^2-sy(inds).^2,z(inds),'smoothingspline','Normalize','on','SmoothingParam',0.95);
indgood=fitpsx(ds)>zrange(1)&fitpsx(ds)<zrange(2);
ds2range=[min(ds(indgood)) max(ds(indgood))];

if nargin>4
    
    plot(ax,z(inds),sx(inds).^2-sy(inds).^2,'.')
    hold(ax,'on')
    sxsort=sort(sx.^2-sy.^2);
    zsort=feval(fitpsx,sxsort);

    plot(ax,zsort,sxsort,'k')
    plot(ax,zrange,[ds2range(1) ds2range(1)],zrange,[ds2range(2) ds2range(2)])
    ylabel(ax,'sx^2-sy^2')
    xlabel(ax,'z (nm)')
     ylim(ax,[q(1) q(2)]);
     title(ax,'Calibration for Gaussian fit')
end

else
    fitpsx=[];
end
out.function=fitpsx;
out.ds2range=ds2range;
end
