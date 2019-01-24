sites=g.locData.SE.sites;

% dt=getFieldAsVectorInd(sites,'evaluation.NPCgeomtryQuantify.templatefit.fitted',4);
% ind = dt>20&dt<100;
% sites = sites(ind);
ff='%2.1f';
f=figure(88);
if isfield(sites(1).evaluation.NPCgeomtryQuantify,'profile') %z-data is there

z0=getFieldAsVector(sites,'evaluation.NPCgeomtryQuantify.profile.Gaussfit.b');
sigma=getFieldAsVector(sites,'evaluation.NPCgeomtryQuantify.profile.Gaussfit.c');
d=getFieldAsVector(sites,'evaluation.NPCgeomtryQuantify.profile.Gaussfit.d');
draw=getFieldAsVector(sites,'evaluation.NPCgeomtryQuantify.profile.Gaussfitraw.d');

subplot(2,4,5)
n=0:1:25;
histogram(sigma,n); xlabel('sigma (nm)')
title(['sigma z: ' num2str(mean(sigma),ff) '\pm' num2str(std(sigma),ff)])


n=20:1:80;
subplot(2,4,6);histogram(abs(d),n); xlabel('d (nm)')
hn=histcounts(abs(d),n);
nf=n(1:end-1)+(n(2)-n(1))/2;
fitp=fit(nf',hn','gauss1');
hold on
plot(nf,fitp(nf))
hold off
title(['d2G: med ' num2str(median(d),ff) ', mean ' num2str(mean(d),ff) '\pm' num2str(std(d),ff) ', fit ' num2str(fitp.b1,ff) '\pm' num2str(fitp.c1/sqrt(2),ff)] )

subplot(2,4,4);histogram(abs(draw),n); xlabel('d (nm)')

hn=histcounts(abs(draw),n);
nf=n(1:end-1)+(n(2)-n(1))/2;
fitp=fit(nf',hn','gauss1');
hold on
plot(nf,fitp(nf))
hold off
title(['draw2G: med ' num2str(median(draw),ff) ', mean ' num2str(mean(draw),ff) '\pm' num2str(std(draw),ff) ', fit ' num2str(fitp.b1,ff) '\pm' num2str(fitp.c1/sqrt(2),ff)] )

dt=getFieldAsVectorInd(sites,'evaluation.NPCgeomtryQuantify.templatefit.fitted',4);
subplot(2,4,7)
histogram(abs(dt),n); xlabel('d (nm)')
title(['template d: median' num2str(median(dt),ff) '\pm' num2str(std(dt),ff)])

hn=histcounts(abs(dt),n);
nf=n(1:end-1)+(n(2)-n(1))/2;
fitp=fit(nf',hn','gauss1');
hold on
plot(nf,fitp(nf))
hold off
title(['dtemplate: med ' num2str(median(dt),ff) ', mean ' num2str(mean(dt),ff) '\pm' num2str(std(dt),ff) ', fit ' num2str(fitp.b1,ff) '\pm' num2str(fitp.c1/sqrt(2),ff)] )

ind = dt<80&dt>20;
zt=getFieldAsVectorInd(sites,'evaluation.NPCgeomtryQuantify.templatefit.fitted',3);
dpl=d';
subplot(2,4,8)
hold off
plot(zt(ind),abs(dpl(ind)),'.')
fline=fit(zt(ind),abs(dpl(ind)),'poly1');
hold on
plot(zt,fline(zt),'r')
xlabel('z');ylabel('distance')
title(['d(z=0) fit: ' num2str(fline.p2,ff) ' Corr:' num2str(corr(abs(d(ind)'),zt(ind')))]);
ylim([25 100])
xlim([-200 200])
end

R0=getFieldAsVector(sites,'evaluation.NPCgeomtryQuantify.Rfit');
subplot(2,4,1)
ff='%2.3f';
% figure(112)
hold off
rn=40:0.5:65;
histogram(abs(R0),rn); xlabel('R (nm)')
title(['fitted radius: ' num2str(mean(R0),ff) '\pm' num2str(std(R0),ff) '\pm' num2str(std(R0)/length(R0),2)])
xlabel('radius (nm)')
xlim([45 60])
ylabel('counts')
hh=histcounts(R0,rn);
fitp=fit(rn(1:end-1)'+(rn(2)-rn(1))/2,hh','gauss1');
hold on
plot(rn,fitp(rn))
fitp

aca=0;aca1=0;aca2=0;acc12=0;
for k=1:length(sites)
    if isfield(sites(k).evaluation.NPCgeomtryQuantify.angular,'actheta')
    ach=sites(k).evaluation.NPCgeomtryQuantify.angular.actheta;
    aca=ach+aca;
    end
    if isfield(sites(k).evaluation.NPCgeomtryQuantify.angular,'ac1')
        aca1=aca1+sites(k).evaluation.NPCgeomtryQuantify.angular.ac1;
        aca2=aca2+sites(k).evaluation.NPCgeomtryQuantify.angular.ac2;
        acc12=acc12+sites(k).evaluation.NPCgeomtryQuantify.angular.cc12;
    end
end
aca=aca/length(sites);
tn=sites(1).evaluation.NPCgeomtryQuantify.angular.thetan;
subplot(2,4,2:3)
aca=aca/2;
% norm=length(tn)-(1:length(tn));
hold off
midp=(length(tn)+1)/2;
dtn=tn(2)-tn(1);
tnn=[tn(midp:end)-2*pi-dtn tn(1:midp-1)];
tnn=tnn/pi*180;
acan=[aca(midp:end) aca(1:midp-1)];
% plot(tn(2:end)*50,aca(2:end)/5);
plot(tnn,acan);
ampm=max(aca(2:end));
ampmin=min(aca);
txtshift=[];
if numel(aca1>0)
    aca1=aca1/length(sites);
    aca2=aca2/length(sites);
    acc12=acc12/length(sites);
    aca2n=[aca2(midp:end) aca2(1:midp-1)];
    aca1n=[aca1(midp:end) aca1(1:midp-1)];
    acc12n=[acc12(midp:end) acc12(1:midp-1)];
    hold on
    plot(tnn,aca1n);
    plot(tnn,aca2n);
    plot(tnn,acc12n);
  
    ampm=max(ampm,max(aca1(2:end)));
    ampm=max(ampm,max(aca2(2:end)));
    ampm=max(ampm,max(acc12(2:end)));
   
    
    %derermine position of maxima of CC
    [~,indm]=max(acc12n);
%     win=10;
%     th=tnn(indm-win:indm+win);cch=acc12n(indm-win:indm+win);
%     fitp=fit(th',cch','poly2');
%     dth=-fitp.p2/fitp.p1/2;
%     plot(th,fitp(th));
    
    win=15;
    th=tnn(indm-win:indm+win);cch=acc12n(indm-win:indm+win);
    ft=fittype('b+a*cos(2*pi*(x-d)/45)');
    fitps=fit(th',cch',ft,'StartPoint',[ampm-ampmin, ampmin,tnn(indm)]);
    plot(th,fitps(th),'r')

    legend('all','ring1','ring2','cross-corr','cos fit')
    win=15;
    indm=indm+45; %next peak
    th=tnn(indm-win:indm+win);cch=acc12n(indm-win:indm+win);
    ft=fittype('b+a*cos(2*pi*(x-d)/45)');
    fitps1=fit(th',cch',ft,'StartPoint',[ampm-ampmin, ampmin,tnn(indm)]);
    indm=round(fitps1.d+length(tnn)/2); %next peak
    th=tnn(indm-win:indm+win);cch=acc12n(indm-win:indm+win);
    fitps1=fit(th',cch',ft,'StartPoint',[ampm-ampmin, ampmin,tnn(indm)]);    
    plot(th,fitps1(th),'r')

    indm=round(fitps1.d+length(tnn)/2)-135; %next peak
    th=tnn(indm-win:indm+win);cch=acc12n(indm-win:indm+win);
    fitps2=fit(th',cch',ft,'StartPoint',[ampm-ampmin, ampmin,tnn(indm)]);
    indm=round(fitps2.d+length(tnn)/2); %next peak
    th=tnn(indm-win:indm+win);cch=acc12n(indm-win:indm+win);
    fitps2=fit(th',cch',ft,'StartPoint',[ampm-ampmin, ampmin,tnn(indm)]);    
    plot(th,fitps2(th),'r')
    
     indm=round(fitps2.d+length(tnn)/2)+45; %next peak
    th=tnn(indm-win:indm+win);cch=acc12n(indm-win:indm+win);
    fitps3=fit(th',cch',ft,'StartPoint',[ampm-ampmin, ampmin,tnn(indm)]);
    indm=round(fitps3.d+length(tnn)/2); %next peak
    th=tnn(indm-win:indm+win);cch=acc12n(indm-win:indm+win);
    fitps3=fit(th',cch',ft,'StartPoint',[ampm-ampmin, ampmin,tnn(indm)]);    
    plot(th,fitps3(th),'r')
    
    ff='%2.1f';
    
     txtshift=[', shift: cos fit1: ' num2str(fitps.d,ff) '°, fit2: ' num2str(fitps1.d-45,ff) '°, fit3: ' num2str(fitps2.d+90,ff) '°, fit4: ' num2str(fitps3.d+45,ff) '°'];
     
     
      plot([0 0],[ampmin ampm*1.1],'k')
end
ylim([ampmin ampm*1.1]);
xlim([tnn(1) tnn(end)])
ax=gca;
% ax.XTick=-180:45:180;
ax.XAxis.MinorTickValues=-180:45/4:180;
ax.XAxis.TickValues=-180:45:180;
ax.XMinorTick='on';
ax.XGrid='on';
xlabel('angle (theta) ')
ylabel('auto/cross correlation')
title(['angular correlation. Nlocs=' num2str(length(d)) txtshift])