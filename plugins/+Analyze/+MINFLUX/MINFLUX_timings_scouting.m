
file=g.locData.loc.filenumber;
tsum=zeros(max(g.locData.loc.itr)+1,1);
maxtime=0;
tval=0;
dtscout=[];
for f=1:max(file)
    indf=file==f;
    time=g.locData.loc.time(indf);
    dt=diff([0;time]);

    tid=g.locData.loc.tid(indf);
    itr=g.locData.loc.itr(indf);
    valid=g.locData.loc.vld(indf);
    
    indval=itr==max(itr)&valid;
    tval=tval+sum(dt(indval));

    for k=0:max(itr)
        tsum(k+1)=tsum(k+1)+sum(dt(itr==k));
    end
    maxtime=maxtime+time(end);
    dtscout=vertcat(dtscout, dt(itr==0));
end
%
figure(88)
tiledlayout(1,2,"TileSpacing","tight")
nexttile

bar(0:max(itr),tsum/maxtime)
xlabel("itr")
ylabel("time (fraction)")
title("itr0: " + num2str(tsum(1)/maxtime*100,3) +"%, valid itr" + max(itr)+": " + num2str(tval/maxtime*100,2) +"%")

nexttile
dtscout=dtscout/1000;
edges = logspace(log10(min(dtscout)), log10(max(dtscout)), 31);
histogram(dtscout, edges)
set(gca, 'XScale', 'log')
xlabel("scouting time (s)")



