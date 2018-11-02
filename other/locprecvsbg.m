 
PSFxpix=1;
EMexcess=1;
PSFypix=1;
photall=[100 200 500 1000 2000 5000 10000];
bg=0:1:1000;
figure(77);
hold off
for k=1:length(photall)
    phot=photall(k);
lpthompson=sqrt((PSFxpix.*PSFypix+1/12)./(phot/EMexcess)+8*pi*(PSFxpix.*PSFypix).^2.*bg./( phot/EMexcess).^2);
           
semilogx(bg,lpthompson/lpthompson(1))
hold on
leg{k}=num2str(phot);
end
legend(leg)
ylim([1 3])
xlabel('background (photons)');
ylabel('localization precision normalized to bg=0')
title('localization precision vs background')