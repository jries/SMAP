sites=g.locData.SE.sites;
z0=getFieldAsVector(sites,'evaluation.NPCgeomtryQuantify.Gaussfit.b');
sigma=getFieldAsVector(sites,'evaluation.NPCgeomtryQuantify.Gaussfit.c');
d=getFieldAsVector(sites,'evaluation.NPCgeomtryQuantify.Gaussfit.d');
f=figure(88);histogram(sigma); xlabel('sigma (nm)')

title(mean(sigma))

figure(89);histogram(d); xlabel('d (nm)')
title(mean(d))