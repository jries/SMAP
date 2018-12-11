sites=g.locData.SE.sites;
z0=getFieldAsVector(sites,'evaluation.NPCgeomtryQuantify.profile.Gaussfit.b');
sigma=getFieldAsVector(sites,'evaluation.NPCgeomtryQuantify.profile.Gaussfit.c');
d=getFieldAsVector(sites,'evaluation.NPCgeomtryQuantify.profile.Gaussfit.d');
ff='%2.1f';
f=figure(88);
subplot(2,4,1)
histogram(sigma); xlabel('sigma (nm)')
title(['profile: sigma z: ' num2str(mean(sigma),ff) '\pm' num2str(std(sigma),ff)])

subplot(2,4,2);histogram(abs(d)); xlabel('d (nm)')
title(['profile: distance: ' num2str(mean(d),ff) '\pm' num2str(std(d),ff)])

dt=getFieldAsVectorInd(sites,'evaluation.NPCgeomtryQuantify.templatefit.fitted',4);
subplot(2,4,3)
histogram(abs(dt)); xlabel('d (nm)')
title(['template d: ' num2str(mean(dt),ff) '\pm' num2str(std(dt),ff)])

R0=getFieldAsVector(sites,'evaluation.NPCgeomtryQuantify.Rfit');
subplot(2,4,4)
histogram(abs(R0)); xlabel('R (nm)')
title(['fitted radius: ' num2str(mean(R0),ff) '\pm' num2str(std(R0),ff)])

aca=0;
for k=1:length(sites)
    if isfield(sites(k).evaluation.NPCgeomtryQuantify.angular,'actheta')
    ach=sites(k).evaluation.NPCgeomtryQuantify.angular.actheta;
    aca=ach+aca;
    end
end
aca=aca/length(sites);
tn=sites(1).evaluation.NPCgeomtryQuantify.angular.thetan;
subplot(2,4,5:6)
% norm=length(tn)-(1:length(tn));
plot(tn(2:end)*50,aca(2:end)/5);
xlabel('distance (theta*r) (nm)')
ylabel('auto correlation')
title('angular correlation')

zt=getFieldAsVectorInd(sites,'evaluation.NPCgeomtryQuantify.templatefit.fitted',3);
subplot(2,4,7)
hold off
plot(zt,abs(dt),'+')
fline=fit(zt,abs(dt),'poly1');
hold on
plot(zt,fline(zt),'r')
xlabel('z');ylabel('distance')