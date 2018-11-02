global se
for k=1:length(se.sites)
    nloc(k)=se.sites(k).evaluation.generalStatistics.Nlayers;
    nlist(k)=se.sites(k).annotation.list1.value;
    nlistgb(k)=se.sites(k).annotation.list2.value==1;

end
figure(87)
hist(nloc,10);
title([mean(nloc),median(nloc)])

ncorn=10-nlist;

figure(88)
indg=nlistgb==true;
hc=hist(ncorn(indg),1:10);
hcf=hc(1:8);
totn=sum(hcf);
hcf=hcf/sum(hcf);
nc=1:8;
fitp=lsqcurvefit(@labelingNPCcorners,.5,nc,hcf,[],[],[],8,4);
plot(nc,hcf*totn,nc,labelingNPCcorners(fitp,nc,8,4)*totn)
title(fitp)