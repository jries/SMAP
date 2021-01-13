if exist('obj','var')
    sites=obj.locData.SE.sites;
else
sites=g.locData.SE.sites;
end


figure(239);
subplot(2,4,[1 2 5 6]);
hold off
frames=(1:1000)';
photall=0*frames;
ft = fittype('a*0+b*exp(-x*c)');
for k=length(sites):-1:1
    frame=double(sites(k).evaluation.beadanalysis.frame);
    [frame,ia]=unique(frame);
    phot=double(sites(k).evaluation.beadanalysis.phot(ia));
%     phot=phot-mean(phot(end-30:end));
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
title(['t1/2=' num2str(1/fita.c,'%5.0f')])


noise=getFieldAsVector(sites,'evaluation.beadanalysis.noise');
subplot(2,4,4)
photrange=0:1000:max(phot0)*1.1;
histogram(phot0,photrange);
% plot(phot0,noise,'+')
title(['phot (t=0), mean: ' num2str(mean(phot0),'%5.0f')]);
% ylabel('noise')
% 
subplot(2,4,7)
histogram(noise,10)
title(['noise mean: ' num2str(mean(noise),'%5.0f')])
% title('noise')
% plot(phot0,decay,'+')
% xlabel('phot (t=0)');
subplot(2,4,8)
histogram(decay,10)
title(['decay, mean: ' num2str(mean(decay),'%5.0f')])
