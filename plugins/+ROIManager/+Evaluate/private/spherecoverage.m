function fitp=spherecoverage(tc,pc,rc,hax)
        tx=-pi/2:pi/128:pi/2;
        hxx=(hist(pc,tx)).^(1/2);
        hm=max(hxx);    
       startp=[hm -pi/2*.9 pi/2*.9 .1 0.1]; 
       lb=[0 -pi -pi/2 0.01 0.01];
       ub=[inf pi/2 pi pi/8 pi/8];
        fitp=lsqcurvefit(@areasegment,startp,tx,(hxx),lb,ub);
        if ~isempty(hax)
          plot(tx,hxx,tx,areasegment(fitp,tx),'Parent',hax)
          hax.NextPlot='add';
          x1=max(fitp(2),-pi/2);
           x2=min(fitp(3),pi/2);
          plot([1 1]*x1,[0 max(hxx)],'Parent',hax)
          plot([1 1]*x2,[0 max(hxx)],'Parent',hax)
          hax.NextPlot='replace';
        end
end


function h=areasegment(fitp,tx)
h=fitp(1)*cos(tx).*(1+erf((tx-fitp(2))/fitp(4))).*(1-erf((tx-fitp(3))/fitp(5)));
h=(h).^(1/2);
end