%filter out single locs
dxs=40;
dzs=100;
minneighbours=5;

%look at neighbours
dx=50;
dz=150;
dn=5;
shiftxy=1500;
minddualc=1; %in units of dn

locsm0=g.locData.getloc({'xnm','ynm','znm'},'layer',2,'Position','all'); 
locso0=g.locData.getloc({'xnm','ynm','znm'},'layer',1,'Position','roi');

locsmr0=g.locData.getloc({'xnm','ynm','znm'},'layer',2,'Position','roi');

nm=countneighbours23D(locsm0,locsm0,dxs,dzs);
nmr=countneighbours23D(locsmr0,locsmr0,dxs,dzs);
no=countneighbours23D(locso0,locso0,dxs,dzs);

locsm.xnm=locsm0.xnm(nm>minneighbours);
locsm.ynm=locsm0.ynm(nm>minneighbours);
locsm.znm=locsm0.znm(nm>minneighbours);

locsmr.xnm=locsmr0.xnm(nmr>minneighbours);
locsmr.ynm=locsmr0.ynm(nmr>minneighbours);
locsmr.znm=locsmr0.znm(nmr>minneighbours);


locso.xnm=locso0.xnm(no>minneighbours);
locso.ynm=locso0.ynm(no>minneighbours);
locso.znm=locso0.znm(no>minneighbours);

neighbours=countneighbours23D(locsm,locso,dx,dz);
clear leg
figure(88);
% plot(xa,ya,'.',xr,yr,'.')
hold off
n=0:dn:quantile(neighbours,.999);
nplot=n(1:end-1)+dn/2;
hompa=histcounts(neighbours(neighbours>0-1),n);
% bar(nplot+dn/2,hompa,'c');
plot(nplot,hompa,'r','LineWidth',2);
leg{1}=['OmpA: ' num2str((1-sum(hompa(1:minddualc))/sum(hompa(1:end)))*100,'%2.0f') '%'];

%shiftx
locso2=locso; locso2.xnm=locso2.xnm+shiftxy;
neighbours2=countneighbours23D(locsm,locso2,dx,dz);
h2=histcounts(neighbours2,n);
hold on 
plot(nplot,h2)
leg{2}=['Shift x: ' num2str((1-sum(h2(1:minddualc))/sum(h2(1:end)))*100,'%2.0f') '%'];

%mirrorx
locso3=locso; locso3.xnm=min(locso.xnm)+max(locso.xnm)-locso.xnm;
neighbours3=countneighbours23D(locsm,locso3,dx,dz);
h3=histcounts(neighbours3,n);
plot(nplot,h3)
leg{3}=['Mirror x: ' num2str((1-sum(h3(1:minddualc))/sum(h3(1:end)))*100,'%2.0f') '%'];

%shifty
locso4=locso; locso4.ynm=locso4.ynm+shiftxy;
neighbours4=countneighbours23D(locsm,locso4,dx,dz);
h4=histcounts(neighbours4,n);
plot(nplot,h4)
leg{4}=['Shift y: ' num2str((1-sum(h4(1:minddualc))/sum(h4(1:end)))*100,'%2.0f') '%'];

%mirrory
locso5=locso; locso5.ynm=min(locso.ynm)+max(locso.ynm)-locso.ynm;
neighbours5=countneighbours23D(locsm,locso5,dx,dz);
h5=histcounts(neighbours5,n);
plot(nplot,h5)
leg{5}=['Mirror y: ' num2str((1-sum(h5(1:minddualc))/sum(h5(1:end)))*100,'%2.0f') '%'];

%mirrorxy
locso6=locso; locso6.ynm=min(locso.ynm)+max(locso.ynm)-locso.ynm;
locso6.xnm=min(locso.xnm)+max(locso.xnm)-locso.xnm;
neighbours6=countneighbours23D(locsm,locso6,dx,dz);
h6=histcounts(neighbours6,n);
plot(nplot,h6)
leg{6}=['Mirror x, y: ' num2str((1-sum(h6(1:minddualc))/sum(h6(1:end)))*100,'%2.0f') '%'];

%shiftmx
locso7=locso; locso7.xnm=locso7.xnm-shiftxy;
neighbours7=countneighbours23D(locsm,locso7,dx,dz);
h7=histcounts(neighbours7,n);
plot(nplot,h7)
leg{7}=['Shift -y: ' num2str((1-sum(h7(1:minddualc))/sum(h7(1:end)))*100,'%2.0f') '%'];
%shifty
locso8=locso; locso8.ynm=locso8.ynm-shiftxy;
neighbours8=countneighbours23D(locsm,locso8,dx,dz);
h8=histcounts(neighbours8,n);
plot(nplot,h8)
leg{8}=['Shift -y: ' num2str((1-sum(h8(1:minddualc))/sum(h8(1:end)))*100,'%2.0f') '%'];

hrand=(h2+h3+h4+h5+h6+h7+h8)/7;
plot(nplot,hrand,'k','LineWidth',2)
leg{9}=['Random: ' num2str((1-sum(hrand(1:minddualc))/sum(hrand(1:end)))*100,'%2.0f') '%'];

ylim([0 hompa(2)*1.1])
xlim([0 nplot(end)])
xlabel('number of Mito neighbours')
ylabel('counts')
legend(leg,'FontSize',22)
ax=gca;
ax.FontSize=18;

[~,fn]=fileparts(g.getPar('lastSMLFile'));
title(fn,'Interpreter','none')
% figure(89)
% plot(locso.xnm,locso.ynm,'.',locso6.xnm,locso6.ynm,'.')


% statistics
ompaM=sum(hompa(minddualc+1:end));
ompaBg=sum(hompa(1:minddualc));
roipos=g.getPar('sr_roiposition');
area=roipos(3)*roipos(4);
locsM=length(locsmr.xnm);
locsO=length(locso.xnm);
ompaMr=round(sum(hrand(minddualc+1:end)));
ompaBgr=round(sum(hrand(1:minddualc)));

sprintf('file \t area (um2) \t ompaMito \t ompaBg \t ompaMitoRand \t ompaBgRand \t locsOmpa \t locsmito')
[~,filen]=fileparts(g.getPar('lastSMLFile'));

outtxt=sprintf([filen '\t' num2str(area) '\t' num2str(ompaM) '\t' num2str(ompaBg) '\t' num2str(ompaMr) '\t' num2str(ompaBgr) '\t' num2str(locsO) '\t' num2str(locsM)]);
clipboard('copy',outtxt)


function neighbours=countneighbours23D(locsm,locso,dx,dz)
sortm=horzcat(locsm.xnm,(1:length(locsm.xnm))');
[sortedm,sortind]=sortrows(sortm);
xr=locsm.xnm(sortind);
yr=locsm.ynm(sortind);
zr=locsm.znm(sortind);

sorto=horzcat(locso.xnm,(1:length(locso.xnm))');
[sortedo,sortindo]=sortrows(sorto);
xa=locso.xnm(sortindo);
ya=locso.ynm(sortindo);
za=locso.znm(sortindo);

                
              

neighbours=countneighbours3Dcirc2(double(xa),double(ya),double(za),double(xr),double(yr),double(zr),double(dx),double(dz));

end