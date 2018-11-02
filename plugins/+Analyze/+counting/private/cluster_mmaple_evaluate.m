function parout=cluster_mmaple_evaluate(par)

if nargin==0
    par.N0_v=10;
    par.pmature_v=0.5;
    par.pblink_v=0.2;
    par.monomer_v=0.2;
    par.ncluster_v=0;
    par.blinkmode.selection='Exponential';
    par.lenbthscale=15;
    par.N0_fit=0;
    par.pmature_fit=1;
    par.pblink_fit=1;
    par.monomer_fit=0;
    par.ncluster_fit=0;
    
end
%evaluate cluster mmaple
global cluster hall hindex



sigmapsf=[cluster(:).meansigmapsf]';
clocs=[cluster(:).locsdf]';
% clocs=[cluster(:).locs]';
ig=clocs>0&sigmapsf>0;
recgui.initaxis(par.resultstabgroup,'scatter')
%     par.resultstabs(3).Title='scatter';
%     ax3=axes('Parent',par.resultstabs(3));
%     axes(ax3)
subplot(2,3,1);
dscatter((([cluster(ig).meansigmapsf]')),[cluster(ig).locsdf]')
xlabel('sigmapsf')
ylabel('locs')

subplot(2,3,2);
dscatter(sqrt(([cluster(ig).stdx].*[cluster(ig).stdy])'),[cluster(ig).locsdf]')
xlabel('stdy')
ylabel('locs')

subplot(2,3,3)
dscatter((([cluster(ig).meanbg]')),[cluster(ig).locsdf]')
xlabel('bg')
ylabel('locs')


subplot(2,3,4)
dscatter((([cluster(ig).meanlocp]')),[cluster(ig).locsdf]')
xlabel('locp')
ylabel('locs')

subplot(2,3,6)
dscatter([cluster(ig).locs]',[cluster(ig).locsdf]')

xlabel('locs')
ylabel('locs df')

% %
% X=[[cluster(ig).locs]',[cluster(ig).meanlocp]',[cluster(ig).meanbg]',[cluster(ig).stdx]',[cluster(ig).stdy]',[cluster(ig).meansigmapsf]'];%,cluster(ig).meanbg',cluster(ig).stdx',cluster(ig).stdy',cluster(ig).meansigmapsf'];
% % [p,score]=pca(X);
% % figure(123)
% % plot(score(:,1),score(:,2),'.')
% idx=kmeans(X,5,'Replicates',10);
% figure(123)
% col=jet(max(idx));
% scatter(X(:,5),X(:,1),12,col(idx,:),'+')


% %filter
minsigma=100;
maxsigma=160;
maxbg=150;

ig=[cluster(:).meansigmapsf]>minsigma & [cluster(:).meansigmapsf]<maxsigma;
% ig=ig&sqrt([cluster(:).stdx].*[cluster(:).stdy])>5 & sqrt([cluster(:).stdx].*[cluster(:).stdy])<35;
% ig=ig&[cluster(:).meanbg]<maxbg;
% ig=ig&[cluster(:).meanlocp]>7 & [cluster(:).meanlocp]<30;

dlocs=[cluster(ig).locsdf]';
% dlocs=[cluster(ig).locs]';
recgui.initaxis(par.resultstabgroup,'histogram')
%     par.resultstabs(4).Title='histogram/fit';
%     ax4=axes('Parent',par.resultstabs(4));
%     axes(ax4)

maxhist=100;
% figure(43)
[y,x]=hist(dlocs,0:maxhist);
sum(y)
length(cluster)
maxfit=20;
minfit=4;
y=y/sum(y(:));
histogram.h=y;
histogram.c=x;
xlim([1,100])




%fit theoretical distribution
% ft = fittype('a*brightnessdistribution( x, 10,1,0 )');
% fit(x,y,ft)
allpar=[par.N0_v,par.pmature_v,par.pblink_v,par.ncluster_v,1];
fitind=[par.N0_fit,par.pmature_fit,par.pblink_fit,par.ncluster_fit,1];

% fitp=lsqcurvefit(@bforfit,startp,x(minfit:maxfit)',y(minfit:maxfit)',[0 0 ],[inf,inf]);
addpar.blinkmode=par.blinkmode.selection;
plot(x,y)
hold on
plot(x,bforfit(allpar,x',addpar),'g')
fitp=fitSome(@bforfit,x(minfit:maxfit)',y(minfit:maxfit)',allpar,fitind,[0 0.01 0.001 0 0],[20,1, 1,2.5, inf],addpar)
% fitp=lsqcurvefit(@bforfit,startp,x(minfit:maxfit)',y(minfit:maxfit)',[0 0.01 0.001 0],[inf,1, 1,2.5]);
plot(x,bforfit(fitp,x',addpar),'r')

hold off
title(fitp);
%assign numlocsreduced to each localization

hall(hindex).histogram=histogram;
hall(hindex).histfit.y=bforfit(fitp,x',addpar);
hall(hindex).fitp=fitp;
% hall(hindex).cluster=cluster;

plotc.x=vertcat(cluster(:).x);
plotc.y=vertcat(cluster(:).y);
plotc.c=vertcat(cluster(:).locsdfv);
plotc.s=vertcat(cluster(:).locprec);

parout.N0_v=10;
parout.pmature_v=fitp(2);
parout.pblink_v=fitp(3);
% parout.monomer_v=0.2;
parout.ncluster_v=fitp(4);

if 0
N0range=6:20;
options= optimoptions('lsqcurvefit');
resnorm=0*N0range;
for k=1:length(N0range)
   startp=[1] ;
   [fitp,resnorm(k)]=lsqcurvefit(@bforfitN0,startp,x(minfit:maxfit)',y(minfit:maxfit)',[0 ],[inf],options,N0range(k));
   fita{k}=fitp;
%    disp([resnorm, N0range(k)])
end


[~,mink]=min(resnorm);
figure(87)
subplot(2,1,1)
plot(N0range,resnorm);
subplot(2,1,2)
plot(x,y)
hold on
plot(x,bforfitN0(fita{mink},x',N0range(mink)),'r')
hold off



% posa.c=posa.numlocsred;
% if 0
    pixrec2=10;
lut=jet;
lut(1,:)=[0 0 0];
pixrec2=5;
imsr=gaussrender(plotc,[0 cluster(1).xmax],[0 cluster(1).ymax],pixrec2, pixrec2,lut,[0 40]);
figure(32)
image(imsr)
colormap jet
colorbar
end
% profile viewer

function fitpout=fitSome(fun,x,y,allpar,fitind,lbi,ubi,addpar)
if nargin<6
    lb=[];
    ub=[];
else
    lb=lbi(fitind==1);
    ub=ubi(fitind==1);
end
startp=allpar(fitind==1);

addp.fitind=fitind==1;
addp.allpar=allpar;
addp.fitfunction=fun;
addp.addpar=addpar;
options= optimoptions('lsqcurvefit');
fitp=lsqcurvefit(@callfitfun,startp,x,y,lb,ub,options,addp);
fitpout=allpar;
fitpout(fitind==1)=fitp;

function y=callfitfun(fitpar,x,addp)
par=addp.allpar;
par(addp.fitind)=fitpar;
addpar=addp.addpar;
y=addp.fitfunction(par,x,addpar);

function y=bforfit(par,x,addpar)

y=par(5)*brightnessdistribution(x,round(par(1)),par(2),par(3),par((4)),'blinkmode',addpar.blinkmode);
% y=par(1)*brightnessdistribution(x,round(par(2)),.50,.33,0,'blinkmode','Exponential');
