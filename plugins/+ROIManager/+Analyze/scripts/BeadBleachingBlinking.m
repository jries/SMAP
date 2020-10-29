sites=g.locData.SE.sites;



figure(239);
subplot(2,4,[1 2 5 6]);
hold off
frames=(1:1000)';
photall=0*frames;
ft = fittype('a*0+b*exp(-x*c)');
for k=length(sites):-1:1
    frame=double(sites(k).evaluation.beadanalysis_2.frame);
    [frame,ia]=unique(frame);
    phot=double(sites(k).evaluation.beadanalysis_2.phot(ia));
    phot=phot-mean(phot(end-30:end));
    plot(frame,phot)
    hold on
    photi=interp1(frame, phot,frames);
    photi(isnan(photi))=phot(1);
    photall=photall+photi;
    phot0(k)=mean(phot(1:10));
    
    startp=[min(phot) max(phot)-min(phot), 1/500];
    fitr=fit(frame,phot,ft,'StartPoint',startp);
    decay(k)=fitr.c;
    
end
ylabel('photons')
xlabel('frames')

subplot(2,4,3)
hold off
plot(frames,photall/k)
startp=[min(photall/k) max(photall/k)-min(photall/k), 1/500];
fita=fit(frames,photall/k,ft,'StartPoint',startp);
hold on
plot(frames,fita(frames));
ylabel('photons')
xlabel('frames')
out.decayall=fita.c;
title(['t1/2=' num2str(1/fita.c)])


noise=getFieldAsVector(sites,'evaluation.beadanalysis_2.noise');
subplot(2,4,4)
photrange=0:1000:max(phot);
histogram(phot,photrange);
% plot(phot0,noise,'+')
title('phot (t=0)');
% ylabel('noise')
% 
subplot(2,4,7)
histogram(noise,10)
title('noise')
% plot(phot0,decay,'+')
% xlabel('phot (t=0)');
subplot(2,4,8)
histogram(decay,10)
title('decay')