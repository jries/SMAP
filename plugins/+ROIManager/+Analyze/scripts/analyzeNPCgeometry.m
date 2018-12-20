sites=g.locData.SE.sites;

% dt=getFieldAsVectorInd(sites,'evaluation.NPCgeomtryQuantify.templatefit.fitted',4);
% ind = dt>20&dt<100;
% sites = sites(ind);

if isfield(sites(1).evaluation.NPCgeomtryQuantify,'profile') %z-data is there

z0=getFieldAsVector(sites,'evaluation.NPCgeomtryQuantify.profile.Gaussfit.b');
sigma=getFieldAsVector(sites,'evaluation.NPCgeomtryQuantify.profile.Gaussfit.c');
d=getFieldAsVector(sites,'evaluation.NPCgeomtryQuantify.profile.Gaussfit.d');
ff='%2.1f';
f=figure(88);
subplot(2,4,5)
histogram(sigma); xlabel('sigma (nm)')
title(['profile: sigma z: ' num2str(mean(sigma),ff) '\pm' num2str(std(sigma),ff)])

subplot(2,4,6);histogram(abs(d)); xlabel('d (nm)')
title(['profile: distance: median' num2str(median(d),ff) '\pm' num2str(std(d),ff)])

dt=getFieldAsVectorInd(sites,'evaluation.NPCgeomtryQuantify.templatefit.fitted',4);
subplot(2,4,7)
histogram(abs(dt)); xlabel('d (nm)')
title(['template d: median' num2str(median(dt),ff) '\pm' num2str(std(dt),ff)])


ind = d<8000000&d>0;
zt=getFieldAsVectorInd(sites,'evaluation.NPCgeomtryQuantify.templatefit.fitted',3);
subplot(2,4,8)
hold off
plot(zt,abs(dt),'+')
fline=fit(zt,abs(dt),'poly1');
hold on
plot(zt,fline(zt),'r')
xlabel('z');ylabel('distance')
title(['d(z=0) fit: ' num2str(fline.p2,ff) ' Corr:' num2str(corr(abs(d(ind)'),zt(ind')))]);
end

R0=getFieldAsVector(sites,'evaluation.NPCgeomtryQuantify.Rfit');
subplot(2,4,1)
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
subplot(2,4,2:3)
% norm=length(tn)-(1:length(tn));
plot(tn(2:end)*50,aca(2:end)/5);
xlabel('distance (theta*r) (nm)')
ylabel('auto correlation')
title('angular correlation')