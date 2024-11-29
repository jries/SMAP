% script to evaluate MINFLUX dc distances over time
% first run cluster_statistics MINFLUX to determine the dual-color tracks
locs=g.locData.getloc({'xnm','ynm','znm','clusterindex','tid','groupindex','time','frame','thi','c_dcind'},'Position','all','grouping','ungrouped');
dcind=locs.c_dcind;
cinds=unique(dcind);
figure(188);
hold off
for k=1:length(cinds)
    indc=dcind==cinds(k);
    indcol1=locs.thi==0;
    indcol2=locs.thi==1;
    xh1=locs.xnm(indc&indcol1);
    yh1=locs.ynm(indc&indcol1);
    th1=locs.time(indc&indcol1);
    xh2=locs.xnm(indc&indcol2);
    yh2=locs.ynm(indc&indcol2);
    th2=locs.time(indc&indcol2);

    indt2=1;
    dx=0*th1;
    dy=0*th1;
    dt=0*th1;
    for t=1:length(th1)
        while indt2<length(th2) && th2(indt2)<th1(t)
            indt2=indt2+1;
        end
        if indt2>1 && th2(indt2)-th1(t)>th1(t)-th2(indt2-1)
            indt2=indt2-1;
        end
        dx(t)=xh1(t)-xh2(indt2);
        dy(t)=yh1(t)-yh2(indt2);
        dt(t)=th1(t)-th2(indt2);
    end


    plot(th1,dx)
    hold on
end
