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
title(['profile: sigma z: ' num2str(mean(sigma),ff) '\pm' num2str(std(sigma),ff)])


n=20:1:70;
subplot(2,4,6);histogram(abs(d),n); xlabel('d (nm)')
title(['profile: distance: median' num2str(median(d),ff) '\pm' num2str(std(d),ff)])


hn=histcounts(abs(d),n);
nf=n(1:end-1)+(n(2)-n(1))/2;
fitp=fit(nf',hn','gauss1');
hold on
plot(nf,fitp(nf))
hold off

subplot(2,4,4);histogram(abs(draw),n); xlabel('d (nm)')
title(['profile: d raw: median' num2str(median(draw),ff) '\pm' num2str(std(draw),ff)])


hn=histcounts(abs(draw),n);
nf=n(1:end-1)+(n(2)-n(1))/2;
fitp=fit(nf',hn','gauss1');
hold on
plot(nf,fitp(nf))
hold off

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
ylim([25 60])
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

aca=0;
for k=1:length(sites)
    if isfield(sites(k).evaluation.NPCgeomtryQuantify.angular,'actheta')
    ach=sites(k).evaluation.NPCgeomtryQuantify.angular.actheta;
    aca=ach+aca;
    end
end
aca=aca/length(sites);
tn=sites(1).evaluation.NPCgeomtryQuantify.angular.thetan;
subplot(2,4,2:3)
% norm=length(tn)-(1:length(tn));
plot(tn(2:end)*50,aca(2:end)/5);
xlabel('distance (theta*r) (nm)')
ylabel('auto correlation')
title(['angular correlation. Nlocs=' num2str(length(d))])